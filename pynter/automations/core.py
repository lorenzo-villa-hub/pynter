#!/usr/bin/env python

import os
from shutil import copyfile
import argparse as ap

from pynter.slurm.job_settings import JobSettings

class Automation:

    def __init__(self, job_script_filename = None, status_filename='exit_status.txt', path=None):
        """
        Abstract class that contains methods useful to automize calculations

        Parameters
        ----------
        job_script_filename : (str), optional
            Name of script for job submission. The default is set in "config.yml".
        status_filename : (str), optional
            Name for status file written after job exit. The default is "exit_status.txt".
        path : (str), optional
            Path of calculation to automize. The default is None. If None current work dir is used
        """        
        
        self.job_script_filename = job_script_filename if job_script_filename else JobSettings().filename
        self.status_filename = status_filename
        self.path = path if path else os.getcwd()


    def copy_files_to_next_step_dir(self,*args):
        """
        Copy files to directory of next step

        Parameters
        ----------
        *args : (str or tuple)
            File names to copy. If is a string the same name is used for root and destination,
            if is a tuple first string is root, second is destination.
        """        
        next_step_path = self.get_next_step()
        for file in args:
            if isinstance(file,str):
                copyfile(os.path.join(self.path,file),os.path.join(next_step_path,file))
            elif isinstance(file, tuple):
                copyfile(os.path.join(self.path,file[0]),os.path.join(next_step_path,file[1]))
            else:
                raise ValueError('Arguments must be tuples of strings')                
        return
        
    def get_common_path(self):
        """
        Get path that is in common to all steps. This method returns the path removing the last directory. 
        """    
        return os.path.normpath(self.path.replace(os.path.basename(self.path),''))


    def get_next_step(self):
        """
        Get path that  of next step directory. 
        """         
        steps_dirs = self.get_ordered_steps()
        current_index = steps_dirs.index(self.path)
        if current_index+1  < len(steps_dirs):
            next_step_dir = steps_dirs[current_index + 1]
            return next_step_dir
        else:
            print('Dir:"%s" - No next step found' %(self.path))
        

    def get_ordered_steps(self):
        """
        Get list of all the directories present if the first subdirectory and orders them.
        The steps must have the step number as first character

        Returns
        -------
        steps_dirs : (list)
            List of paths of all steps in order.
        """
        common_path = self.get_common_path()
        steps_dirs = [f.path for f in os.scandir(common_path) if f.is_dir()]
        steps_dirs.sort() #might wanna add option to change sorting key. in this way the rule of first number works
        
        return steps_dirs
    
    
    def submit_job(self,job_script_path=None,job_script_filename=None):
        """
        Submits job by launching 'sbatch' command on OS

        Parameters
        ----------
        job_script_path : (str), optional
            Path were job script is located. The default is None. If None the next step path is used.
        job_script_filename : (str), optional
            Filename of job script. The default is None. If None the attribute 'job_script_filename' is used.
        """
        wdir = os.getcwd()
        
        if job_script_path is None:
            job_script_path = self.get_next_step()
        if job_script_filename is None:
            job_script_filename = self.job_script_filename
            
        job_script_path = os.path.normpath(job_script_path)
        os.chdir(job_script_path)        
        os.system('sbatch %s' %job_script_filename)
        os.chdir(wdir)       
        return



class CommandHandler(Automation):
    """
    Subclass of Automation to handle commands from command line
    """

    def parser_common_args(self):
        """
        Add general commands that can be used for every program 

        Returns
        -------
        parser : ArgumentParser object
        """

        parser = ap.ArgumentParser()
        
        parser.add_argument('-j','--job-script',help='Job script filename, default is "job.sh"',required=False,default=self.job_script_filename,type=str,metavar='',dest='job_script_filename')
        parser.add_argument('-s','--status',help='Write exit status to file, default is "exit_status.txt"',required=False,default=self.status_filename,type=str,metavar='',dest='status_filename')
        parser.add_argument('-e','--error-check',action='store_true',help='Perform error checking. Default is False',required=False,default=False,dest='error_check')
        
        return parser
    
    def vasp_args(self):
        """
        Add arguments that can be used for VASP

        Returns
        -------
        args : ArgumentParser object
            ArgumentParser object with parsed arguments. An attribute named with every argument "dest" is created.
        """
        
        parser = self.parser_common_args()
        # part for VASP
        parser.add_argument('-c','--contcar',action='store_true',help='Copy CONTCAR to POSCAR of next step',required=False,default=False,dest='contcar')
        parser.add_argument('-W','--wavecar',action='store_true',help='Copy WAVECAR to next step',required=False,default=False,dest='wavecar')
        parser.add_argument('-C','--chgcar',action='store_true',help='Copy CHGCAR to next step',required=False,default=False,dest='chgcar')
        parser.add_argument('-K','--check-kpoints',action='store_true',help='Copy WAVECAR and POSCAR only if KPOINTS of next step are the same',required=False,default=False,dest='check_kpoints')        
           
        args = parser.parse_args()
        
        args.status_filename = args.status_filename if args.status_filename != 'None' else None # if status filename argument is 'None' string set attribute to None
        
        # update class attributes with passed arguments
        self.job_script_filename = args.job_script_filename
        self.status_filename = args.status_filename
        
        return args

        

                    
                  
                    
                    
