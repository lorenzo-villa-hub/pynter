#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:56:45 2020

@author: villa
"""

import os
import re
import warnings

def job_status(path=None,job_script_filename='job_vasp.sh'):
        
    if path == None:
        path = os.getcwd()
    os.system('squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R" > job_list.txt')

    file = os.path.join(path,job_script_filename)
    if os.path.isfile(file):
        job_vasp = open(file)
       
    for line in job_vasp:
        # emulating 'grep' command
        job_line = re.findall(r'#SBATCH --job-name=', line)                    
        if job_line:
            job_name = line.split('=')[1].strip()
            
    
    with open('job_list.txt') as jobs_file:
        
        count = 0
        status_string = None        

        for line in jobs_file:
            split_line = line.split()
            name = split_line[2]
            status = split_line[4]            
 
            if name == job_name:
                
                if status == 'R':
                    status_string = 'RUNNING'
                    count += 1
                elif status == 'PD':
                    status_string = 'PENDING'

        if count > 1:
            warnings.warn(f': More than one job named "{job_name}" are simultaneously running\n',UserWarning)
        
        if not status_string:    
            status_string = 'NOT IN QUEUE'

    os.remove('job_list.txt') 
	        
    return status_string
