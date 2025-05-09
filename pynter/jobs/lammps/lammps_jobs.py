#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 15:37:41 2025

@author: villa
"""


import os
import os.path as op
import warnings
import json
import pandas as pd

from monty.json import jsanitize, MontyDecoder
from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.lammps.outputs import parse_lammps_log

from pynter.jobs.core import Job
from pynter.hpc.slurm import JobSettings
from pynter.lammps.inputs import LammpsInput
from pynter import LOCAL_DIR, SETTINGS



def _parse_lammps_log(path,**kwargs):
    try:
        log = parse_lammps_log(op.join(path,'log.lammps'),**kwargs)
        return log
    except:
        warnings.warn('Reading of log.lammps in "%s" failed'%path)
        return False

input_default_filename = SETTINGS['lammps']['input_filename']
data_default_filename = SETTINGS['lammps']['data_filename']
log_default_filename = SETTINGS['lammps']['log_filename']  


class LammpsJob(Job):

        
    def __init__(self,
                 path=None,
                 inputs=None,
                 job_settings=None,
                 outputs=None,
                 job_script_filename=None,
                 name=None,
                 input_filename=input_default_filename,
                 data_filename=data_default_filename,
                 log_filename=log_default_filename,
                 hostname=None,
                 localdir=None,
                 remotedir=None):
        """
        Class to control and organize inputs and outputs of a LAMMPS job.

        Parameters
        ----------
        path : (str)
            Path where job is stored. If None the work dir is used.
        inputs : (dict)
            Dictionary with input data.
        job_settings : (dict)
            Dictionary with job settings. Documentation in JobSettings class in hpc.slurm module
        outputs : (dict), optional
            Dictionary with output data.
        job_script_filename : (str)
            Filename of job script. The default is taken from the key 'filename' in the job_settings in the config file.
        name : (str)
            Name of the job. If none the name is searched in the job script.
        """
        super().__init__(path=path,
                         inputs=inputs,
                         job_settings=job_settings,
                         outputs=outputs,
                         job_script_filename=job_script_filename,
                         name=name,
                         hostname=hostname,
                         localdir=localdir,
                         remotedir=remotedir)
        
        self._input_filename = input_filename
        self._data_filename = data_filename
        self._log_filename = log_filename
        
        if 'data' in self.inputs.keys():
            self._initial_structure = self.inputs['data'].structure
            
    
    def as_dict(self):        
        """
        Get lammpsjob as dictionary. The Vasprun ouput object is not exported.
        
        Parameters
        ----------
        get_band_structure : (bool), optional
            Export BandStructure as dict. The default is False.
            
        Returns:
            Json-serializable dict representation of lammpsjob.
        """                
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "path_relative":self.path_relative,
             "inputs": jsanitize(self.inputs),             
             "job_settings": self.job_settings.as_dict(),
             "job_script_filename":self.job_script_filename,
             "name":self.name,
             "hostname":self.hostname,
             "localdir":self.localdir,
             "remotedir":self.remotedir}
        
        d["outputs"] = jsanitize(self.outputs)       
        return d


    def to_json(self,path):
        """
        Save lammpsjob object as json string or file

        Parameters
        ----------
        path : (str), optional
            Path to the destination file.  If None a string is exported.
        get_band_structure : (bool), optional
            Export BandStructure as dict. The default is False.

        Returns
        -------
        d : (str)
            If path is not set a string is returned.

        """
        d = self.as_dict()
        if path:
            with open(path,'w') as file:
                json.dump(d,file)
            return
        else:
            return d.__str__()   

    
    @staticmethod
    def from_dict(d):
        """
        Construct LammpsJob object from python dictionary.
        
        Returns
        -------
        lammpsjob object
        """
        #ensure compatibility with old path format
        if 'path_relative' in d.keys() and d['path_relative']:
            path = LOCAL_DIR + d['path_relative']
        elif 'path' in d.keys():
            path = d['path']
        inputs = {}
        if 'inp' in d['inputs'].keys():
            inputs['inp'] = LammpsInput.from_dict(d['inputs']['inp'])
        if 'data' in d['inputs'].keys():
            inputs['data'] = LammpsData.from_dict(d['inputs']['data'])
        job_settings = JobSettings.from_dict(d['job_settings'])
        job_script_filename = d['job_script_filename']
        name = d['name']
        hostname = d['hostname'] if 'hostname' in d.keys() else None
        localdir = d['localdir'] if 'localdir' in d.keys() else None
        remotedir = d['remotedir'] if 'remotedir' in d.keys() else None
        outputs= MontyDecoder().process_decoded(d['outputs'])
        
        lammpsjob = LammpsJob(path=path,
                          inputs=inputs,
                          job_settings=job_settings,
                          outputs=outputs,
                          job_script_filename=job_script_filename,
                          name=name,
                          hostname=hostname,
                          localdir=localdir,
                          remotedir=remotedir)
        
        return lammpsjob
        
    
    @staticmethod
    def from_directory(
                        path=None,
                        job_script_filename=None,
                        input_filename=input_default_filename,
                        data_filename=data_default_filename,
                        log_filename=log_default_filename,
                        atom_style='atomic',
                        load_outputs=True,
                        hostname=None,
                        localdir=None,
                        remotedir=None,
                        **kwargs):
        """
        Builds lammpsjob object from data stored in a directory. Input files are read using Pymatgen VaspInput class.
        Output files are read usign Pymatgen Vasprun class.
        Job settings are read from the job script file.

        Parameters
        ----------
        path : (str)
            Path were job data is stored. If None the current wdir is used. The default is None.
        job_script_filename : (str), optional
            Filename of job script. The default is set in the config file.
        load_outputs : (bool)
            Load job outputs. The default is True.
        kwargs : (dict)
            Arguments to pass to parse_log_lammps.
        Returns
        -------
        lammpsjob object.
        
        """
        path = path if path else os.getcwd()
        inputs = {}
        inputs['inp'] = LammpsInput.from_file(op.join(path,input_filename))
        if os.path.isfile(op.join(path,data_filename)):
            inputs['data'] = LammpsData.from_file(op.join(path,data_filename),atom_style=atom_style)
        outputs = {}
        if load_outputs:
            if op.isfile(op.join(path,log_filename)):
                outputs['log'] = _parse_lammps_log(path,**kwargs)

        job_script_filename = job_script_filename if job_script_filename else SETTINGS['job_script_filename']
        job_settings = JobSettings.from_bash_file(path,filename=job_script_filename)
        
        return LammpsJob(path=path,
                         inputs=inputs,
                         job_settings=job_settings,
                         outputs=outputs,
                         job_script_filename=job_script_filename,
                         input_filename=input_filename,
                         data_filename=data_filename,
                         log_filename=log_filename,
                         hostname=hostname,
                         localdir=localdir,
                         remotedir=remotedir
                         )


    @staticmethod
    def from_json(path_or_string):
        """
        Build LammpsJob object from json file or string.

        Parameters
        ----------
        path_or_string : (str)
            If an existing path to a file is given the object is constructed reading the json file.
            Otherwise it will be read as a string.

        Returns
        -------
        lammpsjob object.

        """
        if op.isfile(path_or_string):
            with open(path_or_string) as file:
                d = json.load(file)
        else:
            d = json.loads(path_or_string)
        return LammpsJob.from_dict(d)


    @property
    def inp(self): # avoid input name in python
        return self.inputs['inp']
    
    @property
    def data(self):
        if 'data' in self.inputs.keys():
            return self.inputs['data']

    @property
    def input_filename(self):
        return self._input_filename

    @property
    def data_filename(self):
        return self._data_filename

    @property
    def log_filename(self):
        return self._log_filename
    
    
    @property
    def formula(self):
        """Complete formula from initial structure (read with Pymatgen)"""
        if self.initial_structure:
            return self.initial_structure.composition.formula
        else:
            return None
    
    @property
    def initial_structure(self):
        if hasattr(self, '_initial_structure'):
            return self._initial_structure
        elif 'data' in self.inputs.keys():
            self._initial_structure = self.inputs['data'].structure
            return self._initial_structure
        else:
            return None      
    
    
    @property
    def log(self):
        if 'log' in self.outputs.keys():
            return self.outputs['log']
        else:
            if not op.exists(op.join(self.path,self.log_filename)):
                warnings.warn(f'"{self.log_filename}" file is not present in Job directory')
            return None

    
    @property
    def is_converged(self):
        if self.outputs:
            if 'convergence' in self.outputs.keys():
                return self.outputs['convergence']
    
    def get_inputs(self,sync=False,atom_style='atomic'):
        if sync:
            self.sync_from_hpc()
        inputs = {}
        inputs['inp'] = LammpsInput.from_file(op.join(self.path,self.input_filename))
        if os.path.isfile(op.join(self.path,self.data_filename)):
            inputs['data'] = LammpsData.from_file(op.join(self.path,self.data_filename),atom_style=atom_style)
        return inputs
    
                
    def get_outputs(self,sync=False,log_filename=log_default_filename,get_output_properties=True,**kwargs):
        """
        Get outputs dictionary from the data stored in the job directory. "log.lammps" is 
        read with Pymatgen

        Parameters
        ----------
        sync : (float)
            Sync data from hpc. The default is False.
        log_filename: (str)
            Name of the lammps log file.
        get_output_properties : (float)
            Parse output properties from log file.
        **kwargs : (dict)
            Arguments for the Vasprun class in pymatgen.
        """
        if sync:
            self.sync_from_hpc()
        self._log_filename = log_filename
        outputs = {}
        if op.isfile(op.join(self.path,log_filename)):
            outputs['log'] = _parse_lammps_log(self.path,**kwargs)

        self.outputs = outputs
        if get_output_properties:
            self.get_output_properties()
        return

        
        

    def get_output_properties(self): #initial idea, to implement further based on specific needs
        self.outputs['convergence'] = self._get_convergence()
        return 
        
    
    def write_input(self):
        self.job_settings.write_bash_file(path=self.path,filename=self.job_settings['filename'])
        self.inp.write_file(
                            filename=op.join(self.path,self.input_filename),
                            ignore_comments=False,
                            keep_stages=True)
        if self.data:
            self.data.write_file(filename=op.join(self.path,self.data_filename))
    
        

    def _get_convergence(self):
        is_converged = None
        if self.outputs:
            if 'log' in self.outputs.keys():
                if self.log is not None:
                    is_converged = False    
                    if self.log and len(self.log) > 1 and type(self.log[-1]) == pd.DataFrame: # check if last element on the list is a dataframe - maybe revisit this?
                        is_converged = True
        return is_converged
        



