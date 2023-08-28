#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 14:23:05 2023

@author: villa
"""

import os
import warnings

from pynter import SETTINGS
from pynter.tools.utils import grep_list

class Slurm:
    
    def __init__(self,**kwargs):
        """
        Object to handle sbatch commands in Slurm. Behaves like a dictionary.

        Parameters
        ----------
        **kwargs : 
            All possible sbatch arguments.
        """
        
        self._arguments, self._arguments_legend = read_possible_slurm_arguments()
        
        defaults = SETTINGS['job_settings']['slurm']
        
        for key,value in defaults.items():
            self.set_argument(key, value)
        
        for key,value in kwargs.items():
            if key in self.arguments_legend.keys():
                key = self.arguments_legend[key]
            self.set_argument(key,value)
            
        
    def __str__(self):
        lines = self.script_lines()
        string = ''.join(lines)
        return string
    
    def __repr__(self):
        return self.settings.__repr__()

    def __len__(self):
        return len(self.settings)

    def __iter__(self):
        return self.settings.keys().__iter__()
    
    def __getitem__(self,key):
        return self.settings[key]
    
    def __setitem__(self,key,value):
        self.set_argument(key,value)
        return
    
    def __eq__(self,other):
        if isinstance(other,str):
            return self.__str__() == other
        elif isinstance(other,Slurm):
            return self.settings == other.settings
        elif isinstance(other,dict):
            return self.settings == other    
        
    def set_argument(self,key,value):
        if key not in self.arguments:
            warnings.warn(f'"{key}" is not listed as possible slurm argument')
        setattr(self, key, value)
        return
    
    
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
    
    @property
    def settings(self):
        settings = {}
        for key,value in self.__dict__.items():
            if key[0] != '_':
                settings.update({key:value})
        settings = dict(sorted(settings.items(), key=lambda item: item[0]))
        return settings
    
    @staticmethod
    def from_string(input_string):
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
     

    def script_lines(self):
        lines = []
        lines.append('#!/bin/sh\n')
        for key,value in self.settings.items():
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