#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:24:40 2023

@author: villa
"""

import os
import yaml

from pynter import get_config_from_default_file, get_vasp_defaults_from_default_file, get_lammps_defaults_from_default_file

def setup_config(subparsers):
    parser_configure = subparsers.add_parser('configure',help='Set up pynter configuration files')
    
    parser_configure.add_argument('-f','--file',help='Configure file ("config.yml","vasp.yml","lammps.yml")',
                               default=None,metavar='',dest='file_to_configure')
    parser_configure.add_argument('-i','--info',action='store_true',help='Print configuration files paths',
                               default=False,dest='info')
    parser_configure.set_defaults(func=create_config)
    
    return
  

def create_config(args):
    run_config(file_to_configure=args.file_to_configure,info=args.info)
    return
  
def run_config(file_to_configure=None,info=False):
    
    filename_config = 'config_v2.yml'
    filename_vasp = 'vasp.yml'
    filename_lammps = 'lammps.yml'
    filepath = os.path.join(os.getenv("HOME"),'.pynter')
    
    # print configuration files path
    if info:
        print('Path of general configuration file: %s' %(os.path.join(filepath,filename_config)))
        print('Path of vasp configuration file: %s' %(os.path.join(filepath,filename_vasp)))
        print('Path of lammps configuration file: %s' %(os.path.join(filepath,filename_lammps)))

        return

    # Set config.yml
    if file_to_configure in [None, filename_config]:
        print('\nSet up configuration file for pynter...\n')
        print('Path of configuration file is: %s' %filepath)
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        
        filename = filename_config
        
        print('Begin setup\n')
        print('HPC setup:')
        hostname = input('Hostname: ')
        localdir = input('Local calculation directory: ')
        localdir = os.path.abspath(localdir) if localdir else None
        remotedir = input('Remote calculation directory: ')
        remotedir = os.path.normpath(remotedir) if remotedir else None
        
        print('\nDefault settings for job script setup: ')
        project_id = input('Project ID: ')
        project_id = project_id if project_id else None
        
        job_script_filename = input('Job script filename (default job.sh): ')
        
        API_KEY = input('\nMaterials Project API_KEY (default None):')
        API_KEY = API_KEY if API_KEY else None
        
 
        job_script_filename = job_script_filename if job_script_filename else 'job.sh'
        
        config = get_config_from_default_file()
        config.update({'HPC': 
                          {'hostname': hostname,
                          'localdir': localdir,
                          'remotedir': remotedir}})
        config.update({'API_KEY':API_KEY})
        
        config['job_script_filename'] = job_script_filename
        config['sbatch']['account'] = project_id
                
        with open(os.path.join(filepath,filename),'w+') as f:
            yaml.dump(config,f) 
            
        print(f'Configuration file {filename} written in {filepath}\n')    
        
        
    # set vasp.yml
    if file_to_configure in [None, filename_vasp]:
    
        filename = filename_vasp
        vasp_default = get_vasp_defaults_from_default_file()
        
        path_vasp_exe = input('Default VASP executable: ')
        path_vasp_exe = path_vasp_exe or None
        
        vasp_default['job_settings_default']['vasp_exe'] = path_vasp_exe
        
        with open(os.path.join(filepath,filename),'w+') as f:
            yaml.dump(vasp_default,f) 
        
        print(f'Vasp default file {filename} written in {filepath}')  

    # set lammps.yml
    if file_to_configure in [None, filename_lammps]:
    
        filename = filename_lammps
        lammps_default = get_lammps_defaults_from_default_file()
        
        input_filename = input('Default LAMMPS input filename (default:input.in):')
        input_filename = input_filename or 'input.in'

        data_filename = input('Default LAMMPS data filename (default:structure.data):')
        data_filename = data_filename or 'structure.data'

        log_filename = input('Default LAMMPS log filename (default:log.lammps):')
        log_filename = log_filename or 'log.lammps'
        
        lammps_default['input_filename'] = input_filename
        lammps_default['data_filename'] = data_filename
        lammps_default['log_filename'] = log_filename


        path_lammps_exe = input('Default LAMMPS executable: ')
        path_lammps_exe = path_lammps_exe or None
        
        lammps_default['job_settings_default']['lammps_exe'] = path_lammps_exe
        
        with open(os.path.join(filepath,filename),'w+') as f:
            yaml.dump(lammps_default,f) 
        
        print(f'Lammps default file {filename} written in {filepath}')  











