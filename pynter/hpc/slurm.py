#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 14:23:05 2023

@author: villa
"""

import os
import os.path as op
import warnings
import copy
import re

from monty.json import MSONable

from pynter import SETTINGS
from pynter.tools.utils import grep_list


def read_possible_sbatch_arguments():
    """
    Read possible sbatch arguments from txt file and map (when present) their short versions 
    (-<letter>) with their long versions (--<word>).

    Returns
    -------
    arguments : (list)
        List containing all possible arguments.
    arguments_legend : (dict)
        Dictionary containing argument mapping.
    """
    path = os.path.abspath(__file__).strip(os.path.basename(__file__))
    filename = 'sbatch_arguments.txt'
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


class Sbatch(dict,MSONable):
    
    def __init__(self,load_defaults=True,**kwargs):
        """
        Object to handle sbatch commands in sbatch. Behaves like a dictionary.

        Parameters
        ----------
        **kwargs : 
            All possible sbatch arguments.
        """
        super().__init__()
        self._arguments, self._arguments_legend = read_possible_sbatch_arguments()
        
        defaults = SETTINGS['sbatch']
        
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
            warnings.warn(f'"{key}" is not listed as possible sbatch argument')
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
    def from_bash_script(input_string,load_defaults=False):
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
            
        return Sbatch(**kwargs,load_defaults=load_defaults)
     

    def get_bash_script_lines(self):
        """
        Get list of lines containing #SBATCH arguments for bash script
        """
        lines = []
        lines.append('#!/bin/sh\n')
        for key,value in self.items():
            if value is not None:
                printed_value = f'={value}' if value else '' 
                lines.append(f'#SBATCH --{key}{printed_value}')
        return lines
                


class JobSettings(dict,MSONable):
    
    def __init__(self,load_sbatch_defaults=True,sbatch={},filename=None,script_lines=[]):
        """
        Handles settings for creation of bash scripts for Job submissions (sbatch). 
        Behaves like a dictionary. The sbatch keys names are consistent with the 
        arguments of the sbatch command on the slurm system. Additional text in the
        script is included as a list of lines. 

        Parameters
        ----------
        load_sbatch_defaults : (bool)
            Load default sbatch arguments from configuration file.
        sbatch : (dict or Sbatch)
            Dictionary with sbatch commands.
        filename : (str)
            Default name for bash file. If None the default is taken from configuration file.
        script_lines : (list)
            Additional lines for the bash script.
        """
        
        super().__init__()
        self._filename = filename or SETTINGS['job_script_filename']
        self._script_lines = script_lines
        
        if sbatch:
            if type(sbatch) is dict:
                sbatch = Sbatch(load_defaults=load_sbatch_defaults,**sbatch)
            elif sbatch.__class__.__name__ == 'Sbatch':
                sbatch = sbatch
            else:
                raise ValueError('sbatch settings must be provided either as dictionary or as Sbatch object')                
        elif not sbatch:
            sbatch = Sbatch(load_defaults=load_sbatch_defaults,**sbatch)
        
        self.sbatch = sbatch
        settings = {} 
        settings['sbatch'] = sbatch
        settings['filename'] = self._filename 
        settings['script_lines'] = self._script_lines
        self.update(settings)
               

    def __repr__(self):
        return self.get_bash_script()
    
    def __print__(self):
        return self.__repr__()


    def __getitem__(self,key):
        if key not in self.keys():
            return self.sbatch[key]
        else:
            return super().__getitem__(key)
    
    def __setitem__(self,key,value):
        args = ['filename','script_lines'] # temporary
        if key not in args:
            self.sbatch.__setitem__(key, value)
        elif key != 'sbatch':
            key = '_' + key #hidden attributes
            setattr(self,key,value)
            super().__setitem__(key,value)
        return

    @property 
    def filename(self):
        return self._filename
    
    @property
    def script_lines(self):
        return self._script_lines

    
    def as_dict(self):
        d = dict(self)
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
        return d       
    
    def copy(self):        
        return copy.deepcopy(self)
    
    @classmethod
    def from_dict(cls,d):
        return cls(**{k: v for k, v in d.items() if k not in ["@module", "@class"]})


    @staticmethod
    def from_bash_file(path,filename='job.sh'):
        """
        Create JobSettings object from file. Cannot read added lines in header and body
        """
        file = op.join(path,filename)
        with open(file) as f:
            input_string = f.read()
        job_settings = JobSettings.from_bash_script(input_string)
        job_settings['filename'] = filename
        return job_settings

    @staticmethod
    def from_bash_script(input_string):
        """
        Initialize object from a bash script
        """
        sbatch = Sbatch.from_bash_script(input_string)
        lines = input_string.split('\n')
        script_lines = []
        for line in lines:
            if all(target_string not in line for target_string in ['#SBATCH','#!']):
                   if line:
                       script_lines.append(line)
         
        return JobSettings(load_sbatch_defaults=False,sbatch=sbatch,script_lines=script_lines)


    def get_bash_script(self):
        lines = self.sbatch.get_bash_script_lines()
        lines.append('')
        lines += self.script_lines 
        return '\n'.join(lines)
    
    
    def replace_line(self,old_string,new_string):
        for index, line in enumerate(self.script_lines):
            if old_string in line:
                new_line = line.replace(old_string, new_string)
                self._script_lines[index] = new_line
                
        return

    
    def write_bash_file(self,path=None,filename=None):
        """
        Write job script file.
        
        Parameters
        ----------
        path : (str), optional
            Path to write job script to. The default is None. If None work dir is used.
        filename : (str), optional
            Filename. If None self.filename is used. The default is None.
        """
        if path:
            if not os.path.exists(path):
                os.makedirs(path)
        filename = filename if filename else self.filename
        complete_path = os.path.join(path,self.filename) if path else filename      
        with open(complete_path,'w') as f:
            f.write(self.get_bash_script())        
        return


