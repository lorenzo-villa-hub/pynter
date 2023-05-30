#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 19:25:00 2020

@author: villa
"""
import os
import os.path as op
from glob import glob
from shutil import copyfile
from pynter.automations.core import CommandHandler, Automation
from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.outputs import Vasprun, Oszicar, Outcar
from pynter.tools.utils import grep


class Base(Automation):
    """
    Subclass of Automation that contains methods for automations of VASP calculations. Works as base class for schemes.
    """ 

    def __init__(self,job_script_filename=None,status_filename='exit_status.txt',path=None):
        super().__init__(job_script_filename,status_filename,path)


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
            print('"vasprun.xml" could not be read, calculation probably did not converge')
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
        dirs.sort()
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

 
    def read_oszicar(self,path=None):
        """
        Get Pymatgen OSZICAR object by reading "OSZICAR" file        

        Parameters
        ----------
        path : (str), optional
            Path of "OSZICAR" file. The default is None. If None "path" attribute is used.

        Returns
        -------
        Oszicar object
        """
        path = path if path else self.path
        return Oszicar(os.path.join(path,'OSZICAR'))


    def read_outcar(self,path=None):
        """
        Get Pymatgen OSZICAR object by reading "OUTCAR" file        

        Parameters
        ----------
        path : (str), optional
            Path of "OUTCAR" file. The default is None. If None "path" attribute is used.

        Returns
        -------
        Outcar object
        """
        path = path if path else self.path
        return Outcar(os.path.join(path,'OUTCAR'))

    
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





class Schemes(Base):
    
    def __init__(self,job_script_filename=None,status_filename='exit_status.txt',path=None,status=[],**kwargs):
        """
        Class to organize automation schemes for Vasp

        Parameters
        ----------
        vasp_automation : (VaspAutomation object)
            VaspAutomation object.
        status : (list), optional
            List of lines containing the report of the automation. The default is [].
        **kwargs :
            If None the default arguments created by the CommandHandler class are used.
        """
        super().__init__(job_script_filename,status_filename,path)
        
        self.status = status

        if kwargs: # custom args
            for key, value in kwargs.items():
                setattr(self,key,value)
        
        else:
            args = CommandHandler().vasp_args() # default args
            for key, value in args.__dict__.items():
                setattr(self,key,value)

    
    def compare_next_step_kpoints(self):
        """
        Compare KPOINTS file in current and next step directory 

        Returns
        -------
        kpoints_are_same : bool
            True if KPOINTS files are the same.
        """
        kpoints_are_same = self.compare_kpoints()
        if kpoints_are_same:
            self.status.append('KPOINTS of current and next step are the same')
        else:
            self.status.append('KPOINTS of current and next step are different')
        return kpoints_are_same 
            
      
    def check_convergence(self):
        """
        Get information on electronic and ionic convergence and adds it to status.

        Returns
        -------
        conv_el : (bool)
            Electronic convergence reached.
        conv_ionic : (bool)
            Ionic convergence reached.

        """        
        conv_el, conv_ionic = self.convergence()
        self.status.append('Job exit. Analysing "vasprun.xml"... ')
        self.status.append(f'\nElectronic convergence: {conv_el}')
        self.status.append(f'Ionic convergence: {conv_ionic}\n') 
        if conv_el == False and conv_ionic == False:
            self.status.append('Calculation probably did not converge, control manually')
        return conv_el,conv_ionic


    def next_step_relaxation_schemes(self):
        """
        If next step exists, transfer output files from the current calculation to 
        the next step calculation and submits the new job. \n
        If "check_kpoints" attribute is set to True, CHGCAR and WAVECAR are copied 
        only if KPOINTS file of current and next step are the same. \n
        CHGCAR,WAVECAR and CONTCAR can be copied if their related attributes ("chgcar","wavecar","contcar") are set to True. \n
        This can be done either manually or by including the commands (-C or --chgcar, -W or --wavecar, -c or --contcar) in the command line.
        """
        if self.get_next_step():
            # if True copy CHGCAR and WAVECAR only if KPOINTS are the same, else copy them in any case
            if self.check_kpoints:
                if self.compare_next_step_kpoints():
                    if self.chgcar:
                        self.transfer_chgcar()
                    if self.wavecar:
                        self.transfer_wavecar()                       
            else:
                if self.chgcar:
                    self.transfer_chgcar()
                if self.wavecar:
                    self.transfer_wavecar()                
            # if True copy CONTCAR to next step POSCAR
            if self.contcar:
                self.transfer_contcar_to_poscar()
                                                  
            self.submit_job()
            self.status.append('\nNext step calculation submitted in dir "../%s"' %os.path.basename(self.get_next_step()))            
        else:
            self.status.append('\nNo next step found, no other calculation submitted')        


    def print_status(self):
        """
        Print status lines on terminal
        """
        print('\n'.join(self.status))


    def step_limits_reached(self):
        """
        Check if electronic or ionic steps limit has been reached. \n
        Electronic steps limit is taken from the "NELM" flag in the INCAR file. \n
        Ionic steps limit is taken from the "NSW" flag in the INCAR file

        Returns
        -------
        el_limit_reached : (bool)
        
        ionic_limit_reached : (bool)
        """
        if self.limit_electronic_steps_reached():
            el_limit_reached = True
            self.status.append('Maximum limit of electronic steps has been reached')
        else:
            el_limit_reached = False
        if self.limit_ionic_steps_reached():
            ionic_limit_reached = True
            self.status.append('Maximum limit of ionic steps has been reached')
        else:
            ionic_limit_reached = False
        return el_limit_reached, ionic_limit_reached
                    

    def resubmit_if_step_limits_reached(self):
        """
        To use if calculation is not converged. Does not perform actions in case the number of files named "out.*"
        in the working directory is >= 2. The idea is that if the calculation hasn't converged after reaching the limit
        twice, is better to control manually. \n
        
        Check if electronic or ionic steps limit has been reached. \n
        In case electronic steps limit is reached, resubmits calculation without doing anything. \n
        In case ionic steps limit is reached copies CONTCAR to POSCAR and resubmits calculation.
        """
        
        number_of_out_files = len(glob(os.path.join(self.path,'out*')) )
        if number_of_out_files < 2:
            self.status.append('Checking if limit of electronic or ionic steps have been reached...')
            el_limit, ionic_limit = self.step_limits_reached()
            
            if el_limit==True and ionic_limit==False:
                self.submit_job(job_script_path = self.path)
                self.status.append('Resubmitting calculation. Be sure "ISTART" in your INCAR is set to 1 or default.')
            
            if el_limit==False and ionic_limit==True:
                copyfile(os.path.join(self.path,'CONTCAR'),os.path.join(self.path,'POSCAR'))
                self.submit_job(job_script_path = self.path)
                self.status.append('Resubmitting calculation. Be sure "ISTART" in your INCAR is set to 1 or default.')
            
            if el_limit==False and ionic_limit==False:
                self.status.append('Problem is not due to steps limit being reached, controll manually...')
        else:
            self.status.append('Two out.* files are already present, control manually...')


    def transfer_chgcar(self):
        """
        Copy CHGCAR file to next step directory
        """
        self.copy_files_to_next_step_dir('CHGCAR')
        self.status.append('CHGCAR copied in dir "../%s"' %os.path.basename(self.get_next_step()))

           
    def transfer_contcar_to_poscar(self):
        """
        Copy CONTCAR file to POSCAR file in next step directory
        """
        self.copy_files_to_next_step_dir(('CONTCAR','POSCAR'))
        self.status.append('CONTCAR copied in POSCAR in dir "../%s"' %os.path.basename(self.get_next_step()))
        return

    
    def transfer_wavecar(self):
        """
        Copy WAVECAR file to next step directory
        """
        self.copy_files_to_next_step_dir('WAVECAR')
        self.status.append('WAVECAR copied in dir "../%s"' %os.path.basename(self.get_next_step()))        
        return
    
                
    def write_status(self,filename=None,path=None):
        """
        Write status lines to file

        Parameters
        ----------
        filename : (str), optional
            Filename. The default is None. If None "status_filename" attribute is used.
        path : (str), optional
            Path of file to write. The default is None. If None "path" attribute is used.
        """

        filename = filename if filename else self.status_filename
        path = path if path else self.path
        if filename:
            with open(filename,'w') as f:
                f.write('\n'.join(self.status))
                print('Automation report written in "%s" file in directory "%s"' %(filename,path))
        else:
            print('Writing status file is disabled, if needed please set "status_filename" kwarg or "--status" optional argument on command line')
            
            
class NEBSchemes(Schemes):


    def clean_NEB_dirs(self,printout=False):
        """
        Clean image directories for NEB calculation. All files but CHGCAR,WAVECAR,POSCAR,CONTCAR and OUTCAR are removed.
        """
        image_dirs = self.find_NEB_dirs()
        for d in image_dirs:
            files = [w[2] for w in os.walk(d)][0]
            for f in files:
                if f not in ('CHGCAR','WAVECAR','POSCAR','OUTCAR','CONTCAR'):
                    os.remove(op.join(d,f))
                    message = f'Removed {op.join(d,f)} \n'
                    if printout:
                        print(message)
                    self.status.append(message)                            
        return
    

    def check_ionic_relaxation_from_outfile(self):
        """
        Useful for NEB because Pymatgen fails to read vasprun file for NEB calculations.
        This function reads the outfile with highest number in the dir and checks for the 
        string: "reached required accuracy - stopping structural energy minimisation". 
        """
        outfiles = glob(os.path.join(self.path,'out*'))
        if outfiles:
            outfiles.sort()
            outfile = outfiles[-1]
            lines = grep('reached required accuracy - stopping structural energy minimisation',outfile)
            if lines:
                print('"reached required accuracy - stopping structural energy minimisation" found in %s' %outfile)
                self.status.append('"reached required accuracy - stopping structural energy minimisation" found in %s' %outfile)
                return True
        else:
            return False
        

    def check_limit_ionic_steps_oszicar(self,printout=False):
        """
        Function to check if limit of ionic steps is reached from OSZICAR file. Reading of vasprun.xml
        failes for NEB calculations. This function checks if the limit is reached from OSZICAR file
        in all images directories. The reference is the NSW tag in the INCAR file in the parent dir.
        """
        limit_reached = True
        image_dirs = self.find_NEB_dirs()
        for d in image_dirs:
            if d != image_dirs[0] and d != image_dirs[-1]:
                image_name = op.basename(d)
                n_steps = len(self.read_oszicar(path=d).ionic_steps)
                nsw = Incar.from_file(op.join(self.path,'INCAR'))['NSW'] # check NSW from INCAR in parent directory
                if printout:
                    self.status.append(f'{image_name}: number of ionic steps:{n_steps}, INCAR["NSW"]:{nsw}')
                if nsw != n_steps:
                    limit_reached = False
        if limit_reached:
            self.status.append('Limit of ionic steps has been reached')
        return limit_reached
            
            
    def check_preconvergence_images(self):
        """
        Check if all SCF calculations of the images are converged
        """
        if self.is_preconvergence():
            convergence = True
            image_dirs = self.find_NEB_dirs()
            for d in image_dirs:
                conv_el, conv_ionic = self.convergence(path=d)
                if conv_el == False or conv_ionic == False:
                    convergence = False
                    self.status.append(f'convergence in {d}: False')
                else:
                    self.status.append(f'convergence in {d}: True')
        else:
            raise ValueError("Current main folder doesn't contain preconvergence calculations of the images")
        
        return convergence
                

    def copy_images_next_step(self):
        """
        Copy CHGCAR,WAVECAR and CONTCAR in POSCAR of respective images in next step folders
        """
        next_step_path = self.get_next_step()
        image_dirs = self.find_NEB_dirs()
        for d in image_dirs:
            image_name = op.basename(d)
            files = [w[2] for w in os.walk(d)][0]
            for f in files:
                if f in ('CHGCAR','WAVECAR','OUTCAR'):
                    source = op.join(d,f)
                    dest = op.join(next_step_path,image_name,f)
                    copyfile(source,dest)
                    self.status.append(f'{f} copied in {dest}')
                if f == 'CONTCAR':
                    source = op.join(d,'CONTCAR')
                    dest = op.join(next_step_path,image_name,'POSCAR')    
                    copyfile(source,dest)
                    self.status.append(f'{f} copied in {dest}')
        return
        

    def copy_images_next_step_and_submit(self):
        """
        Copy CHGCAR,WAVECAR and CONTCAR in POSCAR of respective images in next step folders 
        and submit calculation in next step.
        """
        self.copy_images_next_step()
        self.submit_job()
        self.status.append('Calculation in dir "%s" submitted' %os.path.basename(self.get_next_step()))
        return


    def is_NEB_job_finished(self):
        """
        Returns True if job is converged or if the limit of ionic steps has been reached
        """
        is_job_finished = False
        
        if self.check_ionic_relaxation_from_outfile():
            is_job_finished = True
        elif self.check_limit_ionic_steps_oszicar():
            is_job_finished = True
        
        if is_job_finished:
            self.status.append('NEB job is finished: converged or limit of ionic steps has been reached')
        return is_job_finished
                
                
    def is_preconvergence(self):
        """
        Check if current main folder contains images folders with preconvergence calculations.
        criteria for selection of preconvergence calculation is the presence of the INCAr file
        in the image folder.
        """
        is_preconvergence = False
        image_dirs = self.find_NEB_dirs()
        for d in image_dirs:
            files = [w[2] for w in os.walk(d)][0]
            if 'INCAR' in files:
                is_preconvergence = True
            elif is_preconvergence:
                print('Warning: Inconsistencies in images directories: INCAR file is present but not in all directories')
        return is_preconvergence
        
        
    
        
        
        
        