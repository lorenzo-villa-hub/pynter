#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 17:58:14 2025

@author: lorenzo
"""
import os
import os.path as op
from shutil import copyfile


def are_kpoints_equal(job1,job2):
    """
    Compare KPOINTS btw 2 VaspJob objects.
    """               
    path1 = op.abspath(job1.path)
    path2 = op.abspath(job2.path)
    kpoints1 = job1.kpoints.as_dict() 
    kpoints2 = job2.kpoints.as_dict()
    relative_path1 = op.relpath(path1)
    relative_path2 = op.relpath(path2)
    if kpoints1 == kpoints2:
        message = f'KPOINTS in {relative_path1} and {relative_path2} are the same'
        are_kpts_equal = True
    else:
        message = f'KPOINTS in {relative_path1} and {relative_path2} are different'
        are_kpts_equal = False
    return are_kpts_equal, message

def copy_file_from_jobs(filename,job1,job2):
    """
    Copy file from path of Job1 to path of Job2.
    """
    source = op.join(op.abspath(job1.path),filename)
    destination = op.join(op.abspath(job2.path),filename)
    copyfile(source,destination)
    return f'{op.relpath(source)} copied to {op.relpath(destination)}'
    
def copy_contcar_to_poscar(job1,job2):
    """
    Copy CONTCAR of Job1 to POSCAR of Job2
    """
    path1 = op.abspath(job1.path)
    path2 = op.abspath(job2.path)
    source = op.join(path1,'CONTCAR')
    destination = op.join(path2, 'POSCAR')
    copyfile(source, destination)
    return f'{op.relpath(source)} copied to {op.relpath(destination)}'
    
def is_limit_electronic_steps_reached(job):
    """
    Check if electronic step limit has been reached
    """
    if job.is_converged_electronic is not None:
        ionic_steps = job.vasprun.ionic_steps
        n_electronic_steps = len(ionic_steps[-1]['electronic_steps'])
        nelm = job.vasprun.parameters['NELM']        
        return True if n_electronic_steps == nelm else False   
    else:
        return False

def is_limit_ionic_steps_reached(job):
    """
    Check if ionic step limit has been reached
    """
    if job.is_converged_ionic is not None:
        n_ionic_steps = len(job.vasprun.ionic_steps)
        nsw = job.vasprun.parameters['NSW']
        return True if nsw == n_ionic_steps else False 
    else:
        return False
    
    
    
class VaspAutomation:
    
    def __init__(self,status_filename='exit_status.txt'):
        """
        Class to handle automation workflows for VASP

        Parameters
        ----------
        status_filename: (str)
            Filename to write automation output. The default is 'exit_status.txt'.
        """
        self.status = []
        self.status_filename = status_filename
    
    
    def check_convergence(self,job):
        """
        Check electronic and ionic convergence and update status
        """
        conv_el, conv_ionic = job.is_converged_electronic, job.is_converged_ionic
        self.status.append('Analysing "vasprun.xml"... ')
        self.status.append(f'\nElectronic convergence: {conv_el}')
        self.status.append(f'Ionic convergence: {conv_ionic}\n') 
        if conv_el == False and conv_ionic == False:
            self.status.append('Calculation probably did not converge, control manually')
        return conv_el,conv_ionic
                        
    
    def check_step_limits_reached(self,job):
        """
        Check if electronic or ionic steps limit has been reached. \n
        Electronic steps limit is taken from the "NELM" flag in the INCAR file. \n
        Ionic steps limit is taken from the "NSW" flag in the INCAR file

        Returns
        -------
        el_limit_reached : (bool)
        
        ionic_limit_reached : (bool)
        """
        limit_el_steps_reached, limit_ionic_steps_reached = False, False
        if is_limit_electronic_steps_reached(job):
            limit_el_steps_reached = True
            self.status.append('Maximum limit of electronic steps has been reached')
        if is_limit_ionic_steps_reached(job):
            limit_ionic_steps_reached = True
            self.status.append('Maximum limit of ionic steps has been reached')         
        return limit_el_steps_reached, limit_ionic_steps_reached        
 
    
    def copy_files_from_job1_to_job2(self,
                                     job1,
                                     job2,
                                     filenames,
                                     check_kpoints=False):
        """
        Copy VASP files from directory of Job1 to directory of Job2

        Parameters
        ----------
        job1 : (VaspJob)
            VaspJob object.
        job2 : (VaspJob)
            VaspJob object.
        filenames : (list)
            List of VASP file names to copy. If "CONTCAR" is present, it will be copied 
            as "POSCAR" for Job2.
        check_kpoints : (bool)
            Only copy CHGCAR and WAVECAR if the KPOINTS file are the same btw the two jobs.
        """
        kpoints_equal = True
        if check_kpoints:
            kpoints_equal, kpoints_message = are_kpoints_equal(job1, job2)
            self.status.append(kpoints_message)
        for file in filenames:
            if file=='CONTCAR':
                message = copy_contcar_to_poscar(job1, job2)
                self.status.append(message)
                
            elif file in ['CHGCAR','WAVECAR']:    
                if kpoints_equal or check_kpoints is False:
                    message = copy_file_from_jobs(file, job1, job2)
                    self.status.append(message)
                else:
                    self.status.append(f'Copying of {file} ignored because KPOINTS do not match ')
            
            else:
                message = copy_file_from_jobs(file, job1, job2)
                self.status.append(message)   
    
        
    def write_status(self,path,filename=None):
        """
        Write status lines to file

        Parameters
        ----------
        path : (str), optional
            Path of file to write.
        filename : (str), optional
            Filename. The default is None. If None "status_filename" attribute is used.
        """

        filename = filename or self.status_filename
        if filename:
            with open(filename,'w') as f:
                f.write('\n'.join(self.status))
                print('Automation report written in "%s" file in directory "%s"' %(filename,path))
        else:
            print('Writing status file is disabled, if needed please set "status_filename" kwarg or "--status" optional argument on command line')
            
            
            
class VaspNEBAutomation(VaspAutomation):

    def check_convergence(self,job):
        """
        Same logic as `is_converged` attribute in VaspNEBJob, but returns messages.
        Returns "True" if the calculation is converged,
        or the ionic step limit has been reached reading from the OSZICAR file.
        "False" if reading failed, and "None" if is not present in the outputs dictionary.
        """
        is_converged = job.is_required_accuracy_reached
        if not is_converged:
            self.status.append('Rerquired accuracy NOT reached')
            if job.is_step_limit_reached:
                is_converged = True
                self.status.append('Ionic steps limit ("NSW") reached')
        else:
            self.status.append('Required accuracy reached')
        if is_converged:
            self.status.append('Convergence or step limit reached')
        else:
            self.status.append('Convergence or step limit NOT reached')
        return is_converged
    

    def copy_images_files(self,source_path,dest_path,filenames=['CHGCAR','CONTCAR','WAVECAR']):
        """
        Copy CHGCAR,WAVECAR and CONTCAR in POSCAR of respective images in next step folders
        """
        for image_path in self.find_NEB_dirs(path=source_path):
            image_name = op.basename(image_path)
            for file in filenames:
                if file == 'CONTCAR':
                    source = op.join(image_path,'CONTCAR')
                    dest = op.join(dest_path,image_name,'POSCAR')    
                    copyfile(source,dest)
                else:
                    source = op.join(image_path,file)
                    dest = op.join(dest_path,image_name,file)
                    copyfile(source,dest)                
                self.status.append(f'{op.relpath(source)} copied in {op.relpath(dest)}')
        return


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


    def is_preconvergence(self):
        """
        Check if current main folder contains images folders with preconvergence calculations.
        criteria for selection of preconvergence calculation is the presence of the INCAR file
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
        
        