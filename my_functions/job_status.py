#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:56:45 2020

@author: villa
"""

def job_status(job_script_filename='job_vasp.sh'):
    
    import os
    import re

    os.system('squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R" > job_list.txt')

    if os.path.isfile(job_script_filename):
        job_vasp = open(job_script_filename)
       
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
            print(f'\nWARNING: More than one job named "{job_name}" are simultaneously running\n')
        
        if status_string:    
            print(f'Job "{job_name}" status is: {status_string}')
        else:
            print(f'no job named "{job_name}" has been found in queue')

    os.remove('job_list.txt') 
	        
    return status_string
