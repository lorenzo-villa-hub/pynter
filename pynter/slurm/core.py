#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 14:23:05 2023

@author: villa
"""

import os

def read_possible_slurm_arguments():
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