#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  9 14:54:13 2025

@author: villa
"""

import warnings
from monty.json import MSONable


class LammpsInput(MSONable):
    
    def __init__(self,lines):
        self.lines = lines
    
    def __repr__(self):
        return self.get_string()
    
    def __print__(self):
        return self.__repr__()
    
    def __eq__(self,other):
        self_lines = [line for line in self.lines if line] # remove empty lines
        other_lines = [line for line in other.lines if line] 
        return self_lines == other_lines
    
    def copy(self):
        return LammpsInput(self.lines.copy())
    
    @staticmethod
    def from_string(string):
        return LammpsInput(string.split('\n'))
    
    @staticmethod
    def from_file(filename):
        with open(filename,'r') as file:
            lines = []
            for line in file.readlines():
                lines.append(line.strip('\n'))
        return LammpsInput(lines)
                    
    
    def get_command(self,command):
        return self.get_command_and_indexes(command)[0]

    
    def get_command_and_indexes(self,command):
        target_lines = []
        indexes = []
        for index,line in enumerate(self.lines):
            if command in line:
                target_lines.append(line)
                indexes.append(index)
        if len(target_lines) == 0:
            warnings.warn('No line containing command found')
        return target_lines, indexes
    
    def get_string(self):
        return '\n'.join(self.lines)
        
        
    def insert_line(self,line,indexes=None,previous_command=None):
        if indexes:
            pass
        elif previous_command:
            _, indexes = self.get_command_and_indexes(previous_command)
            if len(indexes) > 1:
                warnings.warn('More than one target line found, ')
        else:
            raise ValueError('either indexes or previous command line need to be provided')
        for idx in indexes:
            self.lines.insert(idx+1,line)
        return
    
    
    def set_command(self,command,new_line):
        _, indexes = self.get_command_and_indexes(command)
        for idx in indexes:
            self.lines[idx] = new_line
        return 
            
    
    def write_file(self,filename):
        with open(filename,'w') as file:
            file.writelines('\n'.join(self.lines))