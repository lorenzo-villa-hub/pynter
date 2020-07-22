#!/usr/bin/env python

import os
import os.path as op
from shutil import copyfile
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import Vasprun
import argparse as ap
from pynter.slurm.job_script import ScriptHandler

class Automation:

    def __init__(self, job_script_filename = None, status_filename='exit_status.txt', path=None):
        """
        Abstract class that contains methods useful to automize calculations

        Parameters
        ----------
        job_script_filename : (str), optional
            Name of script for job submission. The default is 'job.sh'.
        status_filename : (str), optional
            Name for status file written after job exit. The default is 'exit_status.txt'.
        path : (str), optional
            Path of calculation to automize. The default is None. If None current work dir is used
        """        
        
        self.job_script_filename = job_script_filename if job_script_filename else ScriptHandler().filename
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
        
        parser.add_argument('-j','--job-script',help='Job script filename, default is "job_vasp.sh"',required=False,default=self.job_script_filename,type=str,metavar='',dest='job_script_filename')
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
    


class VaspAutomation(Automation):
    """
    Subclass of AutomationHandler that contains methods for automation of VASP calculations.
    """ 

    @property
    def vasprun(self):
        return self.read_vasprun()
    

    
    def compare_kpoints(self,dir1=None,dir2=None):
        """
        Compare KPOINTS files in 2 different folders

        Parameters
        ----------
        dir1 : (str), optional
            First path. The default is None. If None "path" attribute is used.
        dir2 : (str), optional
            Second path. The default is None. If None the next step path is used (determined with "get_next_step()").

        Returns
        -------
        kpoints_are_same : (bool)
            True if KPOINTS files are the same, False if they're different.
        """        
        dir1 = dir1 if dir1 else self.path
        dir2 = dir2 if dir2 else self.get_next_step()       
        current_kpoints = Kpoints().from_file(os.path.join(dir1,'KPOINTS')).as_dict()        
        if self.get_next_step():
            next_kpoints = Kpoints().from_file(os.path.join(dir2,'KPOINTS')).as_dict()
            kpoints_are_same = True if current_kpoints == next_kpoints else False
        else:
            kpoints_are_same = None            
        return kpoints_are_same
    

    def convergence(self,path=None):
        """
        Check electronic and ionic convergence by reading "vasprun.xml" file with Pymatgen.\n
        If reading of vasprun failes, returns False for electronic and ionic convergence.

        Parameters
        ----------
        path : (str), optional
            Path of "vasprun.xml" file. The default is None. If None "path" attribute is used.        

        Returns
        -------
        converged_electronic : (bool)
            Electronic convergence reached.
        converged_ionic : (bool)
            Ionic convergence reached.
        """                
        path = path if path else self.path
        converged_electronic = False
        converged_ionic = False
        try:
            vasprun = self.read_vasprun(path=path)
            converged_electronic = vasprun.converged_electronic
            converged_ionic = vasprun.converged_ionic
        except:
            print('"vasprun.xml" could not be read, calculation probably failed')
            pass
            
        return converged_electronic, converged_ionic
       

    def find_NEB_dirs(self,path=None):
        """
        Find directories of images for NEB calculations. Directories are selected if all characters in the
        directory name are digits.
        """
        dirs = []
        path = path if path else self.path
        path = op.abspath(path)
        for d in os.walk(path):
            directory = d[0]
            image_name = op.relpath(directory,start=path)
            if all(c.isdigit() for c in list(image_name)): #check if folder is image (all characters in folder rel path need to be numbers)
                dirs.append(directory)
        return dirs


    def limit_electronic_steps_reached(self):
        """
        Check if limit of electronic steps has been reached. This methods compares the "NELM" flag 
        in the INCAR file with the number of electronic steps performed in the calculation

        Returns
        -------
        bool
            True if limit reached.
        """
        steps = self.vasprun.ionic_steps
        n_electronic_steps = len(steps[-1]['electronic_steps'])
        nelm = self.vasprun.parameters['NELM']        
        return True if n_electronic_steps == nelm else False
        

    def limit_ionic_steps_reached(self):
        """
        Check if limit of ionic steps has been reached. This methods compares the "NSW" flag 
        in the INCAR file with the number of ionic steps performed in the calculation

        Returns
        -------
        bool
            True if limit reached.
        """
        n_steps = len(self.vasprun.ionic_steps)
        nsw = self.vasprun.parameters['NSW']                
        return True if n_steps == nsw else False

     
    def read_vasprun(self,path=None):
        """
        Get Pymatgen Vasprun object by reading "vasprun.xml"        

        Parameters
        ----------
        path : (str), optional
            Path of "vasprun.xml" file. The default is None. If None "path" attribute is used.

        Returns
        -------
        Vasprun object
        """
        path = path if path else self.path
        return Vasprun(os.path.join(path,'vasprun.xml'))

        

                    
                  
                    
                    
