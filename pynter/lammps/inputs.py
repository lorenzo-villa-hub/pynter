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
    def from_string(string,remove_empty_lines=False):
        """
        Get LammpsInput object from string
        """
        if remove_empty_lines:
            string = "\n".join(line for line in string.splitlines() if line.strip())
        return LammpsInput(string.split('\n'))
    
    @staticmethod
    def from_file(filename):
        """
        Get LammpsInput object from file
        """
        with open(filename,'r') as file:
            lines = []
            for line in file.readlines():
                lines.append(line.strip('\n'))
        return LammpsInput(lines)
                    
    
    def get_command(self,command):
        """
        Get line containing a specific command
        """
        return self.get_command_and_indexes(command)[0]

    
    def get_command_and_indexes(self,command):
        """
        Get line containing a specific command and relative index
        """
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
        """
        Get LAMMPS input string
        """
        return '\n'.join(self.lines)
        
        
    def insert_line(self,line,indexes=None,previous_command=None):
        """
        Insert line at a specific index or after a line containing a target command.

        Parameters
        ----------
        line : (str)
            Line to insert.
        indexes : (int)
            List index in which to insert the line.
        previous_command : (str)
            Insert the line after a specific command. This arg and the indexes
            arg are mutually exclusive, indexes has the priority.
        """
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
        """
        Change a line containing a specific command.

        Parameters
        ----------
        command : (str)
            Target command to reset..
        new_line : (str)
            New line to insert.
        """
        _, indexes = self.get_command_and_indexes(command)
        for idx in indexes:
            self.lines[idx] = new_line
        return 
            
    
    def write_file(self,filename):
        """
        Write LAMMPS input file
        """
        with open(filename,'w') as file:
            file.writelines('\n'.join(self.lines))