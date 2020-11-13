#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 14:04:08 2020

@author: villa
"""
import os
import yaml

print('Set up configuration file for pynter...\n\n')

filepath = os.path.join(os.getenv("HOME"),'.pynter')
print('Path of configuration file is: %s' %filepath)
if not os.path.exists(filepath):
    os.makedirs(filepath)

filename = input('Choose file name (default: config.yml): ')
filename = filename if filename else 'config.yml' 

print('Begin setup\n')
print('HPC setup:')
hostname = input('Hostname: ')
localdir = input('Local calculation directory: ')
workdir = input('Remote calculation directory: ')

print('\nDefault settings for job script setup: ')
project_id = input('Project ID: ')
email = input('email address: ')
path_exe = input('Path of default executable: ')
job_script_filename = input('Job script filename (default job.sh): ')

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

config = {
    
    'HPC': 
          {'hostname': hostname,
          'localdir': localdir,
          'workdir': workdir},
    'API_KEY': 'DSR45TfHVuyuB1WvP1',
    'job_settings': 
          {'project_id': project_id,
          'name': 'no_name',
          'array_size': None,
          'email': email,
          'nodes': 4,
          'cores_per_node': 24,
          'output_filename': 'out.%j',
          'error_filename': 'err.%j',
          'timelimit': '24:00:00',
          'memory_per_cpu': 2400,
          'processor': 'avx2',
          'modules': ['intel/2019.2', 'intel/2019.3', 'intelmpi/2019.3', 'fftw/3.3.8'],
          'path_exe': path_exe,
          'add_stop_array': True,
          'add_automation': None,
          'add_lines_header': None,
          'add_lines_body': None,
          'filename': job_script_filename},
    'dbconfig': 
          {'vasp': 
               {'collection': collection,
               'database': database,
               'host': host,
               'password': None,
               'port': port,
               'user': None}}}

    
with open(os.path.join(filepath,filename),'w+') as f:
    yaml.dump(config,f)
    
    
print(f'Configuration file {filename} written in {filepath}')    
    
    