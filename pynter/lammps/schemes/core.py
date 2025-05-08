#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  8 14:44:42 2025

@author: villa
"""
import os
import os.path as op
import warnings

from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile

from pynter import SETTINGS
from pynter.hpc.slurm import JobSettings
from pynter.jobs.lammps.lammps_jobs import LammpsJob


class DefaultInputs:
    
    def __init__(self,structure=None,sort_structure=True):
        """
        Parameters
        ----------
        structure : Pymatgen Structure object.
        sort_structure (bool) : Sort structure.
        """
        if structure and sort_structure:
            structure.sort()
        self._structure = structure if structure else None 
        defaults = SETTINGS['lammps']
        
        self.default_potcar_symbols =  defaults['default_potcar_symbols']
        if self.structure:
            self.potcar_symbols = [self.default_potcar_symbols[el.symbol] for el in self.structure.composition.elements]
        self.incar_default_flags = defaults['incar_default']


    @property
    def structure(self):
        return self._structure
    
    @structure.setter
    def structure(self,structure):
        self._structure = structure
        
        
    def get_lammps_input_default(input_string=None) 



class DefaultJobSettings:
    
    def __init__(self,
                 lammps_exe=None,
                 lammps_sbatch=None,
                 modules=None,
                 lines_before_srun=None,
                 lines_after_srun=None):
        """
        Class to handle LAMMPS specific part of sbatch submission script.
        Args that are not provided are taken from lammps.yml file in ~/.pynter.

        Parameters
        ----------
        lammps_exe : (str)
            Path to lammps executable. 
        lammps_sbatch : (dict)
            Sbatch dictionary tailored for lammps.
        modules : (list)
            Name of modules to load with "module load".
        lines_before_srun : (list)
            Script lines to be placed before the "srun" command.
        lines_after_srun : (list)
            Script lines to be placed before the "srun" command.
        """
        
        default_settings = SETTINGS['lammps']['job_settings_default']
        self.lammps_exe = lammps_exe or default_settings['lammps_exe']
        self.lammps_sbatch = lammps_sbatch or default_settings['lammps_sbatch']
        self.modules = modules or default_settings['modules']
        self.lines_before_srun = lines_before_srun or default_settings['lines_before_srun']
        self.lines_after_srun = lines_after_srun or default_settings['lines_after_srun']
    
    @property
    def lammps_script_lines(self):
        if self.modules:
            lines = self.get_modules_lines()
        else:
            lines = []
        if self.lines_before_srun:
            lines += self.lines_before_srun
        if self.lammps_exe:
            lines += [self.get_srun_line()]
        if self.lines_after_srun:
            lines += self.lines_after_srun
        return lines
        
    def get_modules_lines(self):
        if self.modules:
            return ['module purge'] + [f'module load {module}' for module in self.modules] 

    def get_srun_line(self):
        return 'srun %s' %self.lammps_exe

    def get_updated_job_settings(self,job_settings):
        """
        Get a copy of a JobSettings object updated with lammps specific settings
        """
        job_settings = job_settings.copy()
        for k,v in self.lammps_sbatch.items():
            job_settings[k] = v
        lines = self.lammps_script_lines
        job_settings['script_lines'] += lines
        return job_settings
    
    

class InputSets:
    
    def __init__(self,path=None,structure=None,input_string=None,
                 job_settings=None,name=None,add_parent_folder=False):
        """
        Parameters
        ----------
        path : (str)
            Main path for the scheme. If None the current work dir is used. The default is None.
        vaspinput : (Pymatgen VaspInput object), optional
            Set of VASP inputs, the default is None. if provided the other inputs of the class (structure,incar_settings,kpoints,potcar) are not needed 
        structure : (Pymatgen Structure object), optional
            Pymatgen Structure object.
        incar_settings : (Dict), optional
            Dictionary with incar flags. The default is None. If None the default settings for PBE functional from the DefaultInputs class are used.
        kpoints : (Pymatgen Kpoints object), optional
            Pymatgen Kpoints object. The default is None. If None the default settings from the DefaultInputs class are used.
        potcar : (Pymatgen kpoints object), optional
            Pymatgen kpoints object. The default is None. If None the default settings from the DefaultInputs class are used.
        job_settings : (Dict), optional
            Dictionary with job settings to create job script, parameters are defined in JobSettings class function. The default is None.\n
            If job_settings is None, the 'job-name' key will be added, the value is the 'name' argument if provided, if 'name' arg is \n
            None the value will be: 'no_name'.
        name : (str), optional
            Name for the system to set up scheme for. The default is None.
        add_parent_folder : (bool), optional
            Add folder to the path names like the name of the InputSets. Default is False.
        """
        self.path = op.abspath(path) if path else os.getcwd()
        if structure:
            self.structure = structure          
        else:
            warnings.warn('Structure is not provided explicitly, make sure is defined in input string')
            
        if job_settings:
            self.job_settings = JobSettings(**job_settings)
        elif name:
            self.job_settings = JobSettings(sbatch={'job-name':name})
        else:
            self.job_settings = JobSettings(sbatch={'job-name':'no_name'})
        self.name = name if name != None else self.job_settings['job-name'] 
        
        if self.name: # this return false if name is '', the previuos line considers only if name is None
            self.job_settings['job-name'] = self.name
        
        self.job_settings = DefaultJobSettings().get_updated_job_settings(self.job_settings)    

        if add_parent_folder:
            self.path = op.join(self.path,self.name)
            
    
    def get_lammpsjob(self,atom_style='atomic',add_to_job_name=None,add_to_job_path=None):
        """
        Generate LammpsJob object from the input settings of the class.

        Parameters
        ----------
        add_to_job_name : (str), optional
            String to be added to 'name' key in job_settings dictionary and to name attribute of LammpsJob.
        add_to_job_path : (str), optional
            String to be added to self.path. The complete path will be input in 'path' arg of LammpsJob.

        Returns
        -------
        LammpsJob object
        """        
        job_settings = self.job_settings.copy()
        
        if add_to_job_name:
            job_settings['job-name'] = '_'.join([job_settings['job-name'],add_to_job_name])        
            jobname = '_'.join([self.name,add_to_job_name])
        else:
            jobname = self.name
            
        if add_to_job_path:
            jobpath = op.join(self.path,add_to_job_path)
        else:
            jobpath = self.path
        
        lammps_input = LammpsInputFile(self.input_string)
        inputs = {'inp':lammps_input}
        if self.structure:
            lammps_data = LammpsData.from_structure(structure=self.structure,atom_style=atom_style)
            inputs['data'] = lammps_data
            
        job_script_filename = job_settings['filename'] if 'filename' in job_settings.keys() else None
        lammpsjob = LammpsJob(path=jobpath,inputs=inputs,job_settings=job_settings,job_script_filename=job_script_filename,name=jobname)
        
        return lammpsjob  