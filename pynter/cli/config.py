#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:24:40 2023

@author: villa
"""

import os
import yaml

from pynter import get_config_from_default_file, get_vasp_defaults_from_default_file

def setup_config(subparsers):
    parser_configure = subparsers.add_parser('configure',help='Set up pynter configuration files')
    
    parser_configure.add_argument('-e','--exclude',help='Exclude file ("config.yml" or "vasp.yml")',
                               default=None,metavar='',dest='exclude')
    parser_configure.add_argument('-i','--info',action='store_true',help='Print configuration files paths',
                               default=False,dest='info')
    parser_configure.set_defaults(func=create_config)
    
    return
  

def create_config(args):
    run_config(exclude=args.exclude,info=args.info)
    return
  
def run_config(exclude=None,info=False):
    
    filename_config = 'config.yml'
    filename_vasp = 'vasp.yml'
    filepath = os.path.join(os.getenv("HOME"),'.pynter')
    
    # print configuration files path
    if info:
        print('Path of general configuration file: %s' %(os.path.join(filepath,filename_config)))
        print('Path of vasp configuration file: %s' %(os.path.join(filepath,filename_vasp)))
        return

    # Set config.yml
    if exclude != filename_config:
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
        workdir = input('Remote calculation directory: ')
        workdir = os.path.normpath(workdir) if workdir else None
        
        print('\nDefault settings for job script setup: ')
        project_id = input('Project ID: ')
        project_id = project_id if project_id else None
        email = input('email address (Default None): ')
        email = email if email else None
        path_exe = input('Path of default executable: ')
        path_exe = path_exe if path_exe else None
        job_script_filename = input('Job script filename (default job.sh): ')
        
        API_KEY = input('\nMaterials Project API_KEY (default None):')
        API_KEY = API_KEY if API_KEY else None
        
        print('\nDatabase configuration:')
        host = input('Host (default: 127.0.0.1 ): ')
        port = input('Port (default: 27017): ')
        database = input('Database (default: pynter): ')
        collection = input('Collection (default vasp): ')
        
        
        job_script_filename = job_script_filename if job_script_filename else 'job.sh'
        host = host if host else '127.0.0.1'
        port = port if port else 27017
        database = database if database else 'pynter'
        collection = collection if collection else 'vasp'
        
        config = get_config_from_default_file()
        config.update({'HPC': 
                          {'hostname': hostname,
                          'localdir': localdir,
                          'workdir': workdir}})
        config.update({'API_KEY':API_KEY})
        
        config['job_settings']['filename'] = job_script_filename
        config['job_settings']['path_exe'] = path_exe
        config['job_settings']['slurm']['account'] = project_id
        config['job_settings']['slurm']['mail-user'] = email
        
        config.update({'dbconfig': 
                          {'vasp': 
                                {'collection': collection,
                                'database': database,
                                'host': host,
                                'password': None,
                                'port': port,
                                'user': None}}})
                
        with open(os.path.join(filepath,filename),'w+') as f:
            yaml.dump(config,f) 
            
        print(f'Configuration file {filename} written in {filepath}\n')    
        
        
    # set vasp.yml
    if exclude != filename_vasp:
    
        filename = filename_vasp
        vasp_default = get_vasp_defaults_from_default_file()
        
        with open(os.path.join(filepath,filename),'w+') as f:
            yaml.dump(vasp_default,f) 
        
        print(f'Vasp default file {filename} written in {filepath}')  













