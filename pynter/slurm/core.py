#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 14:23:05 2023

@author: villa
"""

import os
import warnings
import copy

from monty.json import MSONable

from pynter import SETTINGS
from pynter.tools.utils import grep_list


class Slurm(dict,MSONable):
    
    def __init__(self,load_defaults=True,**kwargs):
        """
        Object to handle sbatch commands in Slurm. Behaves like a dictionary.

        Parameters
        ----------
        **kwargs : 
            All possible sbatch arguments.
        """
        super().__init__()
        self._arguments, self._arguments_legend = read_possible_slurm_arguments()
        
        defaults = SETTINGS['job_settings']['slurm']
        
        settings = {}
        if load_defaults:
            for key,value in defaults.items():
                settings[key] = value
        
        for key,value in kwargs.items():
            if key in self.arguments_legend.keys():
                key = self.arguments_legend[key]
            settings[key] = value
        
        settings = dict(sorted(settings.items(), key=lambda item: item[0]))
        self.update(settings)
        
    
    def __setitem__(self,key,value):
        if key not in self.arguments:
            warnings.warn(f'"{key}" is not listed as possible slurm argument')
        else:
            super().__setitem__(key,value)
        return
    
    def as_dict(self):
        d = dict(self)
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
        return d        
    
    def copy(self):
        return copy.deepcopy(self)
    
    def update(self, other=None, **kwargs):
        if other is not None:
            if hasattr(other, "items"):
                for key,value in other.items():
                    self.__setitem__(key, value)
            else:
                for key, value in other:
                    self.__setitem__(key, value)
        return
    
    @classmethod
    def from_dict(cls,d):
        return cls(**{k: v for k, v in d.items() if k not in ["@module", "@class"]})
    
    
    @property
    def arguments(self):
        """
        List of all possible sbatch arguments
        """
        return self._arguments
    
    @property
    def arguments_legend(self):
        """
        Dictionary which maps the short version (-<letter>) of the
        arguments with the long one (--<word>)
        """
        return self._arguments_legend
    
    
    @staticmethod
    def from_bash_script(input_string):
        """
        Initialize object from a bash script
        """
        lines = input_string.split('\n')
        kwargs = {}
        string = '#SBATCH'
        target_lines = grep_list(string,lines)
        for line in target_lines:  
            words = line.split(' ')
            if '-' in line and '--' not in line:
                if len(words) > 2:
                    pair = words[1].strip('-') , words[2]
                else:
                    pair = words[1].strip('-')
            elif '--' in line:
                pair =  words[1].strip('--').split('=')
            key = pair[0]
            value = pair[1] if len(pair) > 1 else ''
            if value.isdigit():
                value = int(value)
            kwargs.update({key:value})
            
        return Slurm(**kwargs)
     

    def get_bash_script_lines(self):
        """
        Get list of lines containing #SBATCH arguments for bash script
        """
        lines = []
        lines.append('#!/bin/sh\n')
        for key,value in self.items():
            if value is not None:
                printed_value = f'={value}' if value else '' 
                lines.append(f'#SBATCH --{key}{printed_value}\n')
        return lines
                
                
                

def read_possible_slurm_arguments():
    """
    Read possible slurm arguments from txt file and map (when present) their short versions 
    (-<letter>) with their long versions (--<word>).

    Returns
    -------
    arguments : (list)
        List containing all possible arguments.
    arguments_legend : (dict)
        Dictionary containing argument mapping.
    """
    path = os.path.abspath(__file__).strip(os.path.basename(__file__))
    filename = 'slurm_arguments.txt'
    with open(os.path.join(path,filename),'r') as file:
        lines = file.readlines()
  
    arguments = []
    arguments_legend = {}
    for line in lines:
        line = line.strip('\n')
        elements = line.split(',')
        arg = elements[0]
        arguments.append(arg)
        if len(elements) > 1:
            arg_short = elements[1]
            arguments_legend.update({arg_short:arg})

    return arguments, arguments_legend