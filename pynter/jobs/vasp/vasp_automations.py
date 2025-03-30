#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 17:58:14 2025

@author: lorenzo
"""
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
    message = f'KPOINTS in {relative_path1} and {relative_path2} are the same'
    return kpoints1 == kpoints2, message

def copy_file_from_jobs(filename,job1,job2):
    """
    Copy file from path of Job1 to path of Job2.
    """
    source = op.join(op.abspath(job1.path1),filename)
    destination = op.join(op.abspath(job2.path2),filename)
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
        if is_limit_electronic_steps_reached:
            limit_el_steps_reached = True
            self.status.append('Maximum limit of electronic steps has been reached')
        if self.limit_ionic_steps_reached():
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
        for file in filenames:
            if file=='CONTCAR':
                message = copy_contcar_to_poscar(job1, job2)
                self.status.append(message)
                
            elif file in ['CHGCAR','WAVECAR']:
                kpoints_equal, message = are_kpoints_equal(job1, job2)
                self.status.append(message)
                if kpoints_equal or check_kpoints is False:
                    message = copy_file_from_jobs(file, job1, job2)
                    self.status.append(message)
            
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
    
    def copy_images_files_from_job1_to_job2(self,job1,job2,filenames=['CHGCAR','CONTCAR','WAVECAR']):
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
        """
        path2 = job2.path
        for image_path in job1.image_dirs:
            image_name = op.path.basename(image_path)
            for file in filenames:
                if file == 'CONTCAR':
                    source = op.join(image_path,'CONTCAR')
                    dest = op.join(path2,image_name,'POSCAR')    
                    copyfile(source,dest)
                else:
                    source = op.join(image_path,file)
                    dest = op.join(path2,image_name,file)
                    copyfile(source,dest)                
                self.status.append(f'{op.relpath(source)} copied in {op.relpath(dest)}')
        
        
        