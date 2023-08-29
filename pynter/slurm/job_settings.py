#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 15:16:02 2023

@author: villa
"""

from monty.json import MSONable

from pynter import SETTINGS
from pynter.slurm.core import Slurm


class JobSettings(dict,MSONable):
    
    # switch to key value format for sbatch args
    def __init__(self,slurm=None,filename='job.sh',array_size=None,modules=None,path_exe=None,
                 add_stop_array=False,add_automation=False,add_lines_header=None,
                 add_lines_body=None,**kwargs):
        """
        Parameters
        ----------
        slurm : (Slurm) Slurm object. Handles all argments related to #SBATCH. 
            If None the defualt values are used, which will be updated with the
            user-defined **kwargs.
        array_size: (int) Number of jobs for array \n
        modules: (list) List of modules to be loaded
        path_exe: (str) Path to executable \n             
        add_stop_array : (Bool), Add lines for stopping array jobs when calculation is converged. \n                
            If True is effective only if key 'array_size' is not None
        add_automation : (str) , Automation script to add to the file.                
        add_lines_header : (List) , Lines to add in the header part of the file.
        add_lines_body : (List) , Lines to add in the body part of the file.
        **kwargs : dict. Additional kwargs which will be automatically assigned to Slurm.
        """
        super().__init__()
        
        settings = {}
        default_settings = SETTINGS['job_settings']
        
        if slurm:
            if type(slurm) is dict:
                self.slurm = Slurm(**slurm)
            elif slurm.__class__.__name__ != 'Slurm':
                raise ValueError('Slurm settings must be provided either as dictionary or as Slurm object')

        elif not slurm and 'slurm' not in kwargs.keys():
            slurm = Slurm(**kwargs)
        
        settings['slurm'] = slurm
    
        settings['filename'] = filename if filename is not None else default_settings['filename']
        settings['array_size'] = array_size if array_size is not None else default_settings['array_size']
        settings['modules'] = modules if modules is not None else default_settings['modules']
        settings['path_exe'] = path_exe if path_exe is not None else default_settings['path_exe']
        settings['add_stop_array'] = add_stop_array if add_stop_array else default_settings['add_stop_array']
        settings['add_automation'] = add_automation if add_automation else default_settings['add_automation']
        settings['add_lines_header'] = add_lines_header if add_lines_header is not None else default_settings['add_lines_header']
        settings['add_lines_body'] = add_lines_body if add_lines_body is not None else default_settings['add_lines_body']

        for key,value in settings.items():
            setattr(self,key,value)
            
        self.update(settings)
    
    def __getitem__(self,key):
        if key not in self.keys():
            return self.slurm[key]
        else:
            return super().__getitem__(key)
    
    def __setitem__(self,key,value):
        if key not in self.keys():
            self.slurm.__setitem__(key, value)
        else:
            super().__setitem__(key,value)
        return
    
        