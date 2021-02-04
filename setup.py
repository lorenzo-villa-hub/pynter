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
email = input('email address (Default None): ')
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
          'modules': None,
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
    
print(f'Configuration file {filename} written in {filepath}\n')    
    

vasp_default = {
  'incar_default': {'ALGO': 'Normal',
  'AMIX': 0.2,
  'EDIFF': 1e-06,
  'EDIFFG': -0.05,
  'ENCUT': 550,
  'IBRION': 2,
  'ISIF': 2,
  'ISMEAR': 0,
  'ISPIN': 1,
  'LCHARG': '.TRUE.',
  'LORBIT': 10,
  'LREAL': '.FALSE.',
  'LWAVE': '.TRUE.',
  'NELM': 200,
  'NSW': 0,
  'SIGMA': 0.05},
 'default_potcar_symbols': {'Ac': 'Ac',
  'Ag': 'Ag',
  'Al': 'Al',
  'Ar': 'Ar',
  'As': 'As',
  'Au': 'Au',
  'B': 'B',
  'Ba': 'Ba_sv',
  'Be': 'Be_sv',
  'Bi': 'Bi',
  'Br': 'Br',
  'C': 'C',
  'Ca': 'Ca_sv',
  'Cd': 'Cd',
  'Ce': 'Ce',
  'Cl': 'Cl',
  'Co': 'Co',
  'Cr': 'Cr_pv',
  'Cs': 'Cs_sv',
  'Cu': 'Cu_pv',
  'Dy': 'Dy_3',
  'Er': 'Er_3',
  'Eu': 'Eu',
  'F': 'F',
  'Fe': 'Fe_pv',
  'Ga': 'Ga_d',
  'Gd': 'Gd',
  'Ge': 'Ge_d',
  'H': 'H',
  'He': 'He',
  'Hf': 'Hf_pv',
  'Hg': 'Hg',
  'Ho': 'Ho_3',
  'I': 'I',
  'In': 'In_d',
  'Ir': 'Ir',
  'K': 'K_sv',
  'Kr': 'Kr',
  'La': 'La',
  'Li': 'Li_sv',
  'Lu': 'Lu_3',
  'Mg': 'Mg_pv',
  'Mn': 'Mn_pv',
  'Mo': 'Mo_pv',
  'N': 'N',
  'Na': 'Na',
  'Nb': 'Nb_pv',
  'Nd': 'Nd_3',
  'Ne': 'Ne',
  'Ni': 'Ni_pv',
  'Np': 'Np',
  'O': 'O',
  'Os': 'Os_pv',
  'P': 'P',
  'Pa': 'Pa',
  'Pb': 'Pb_d',
  'Pd': 'Pd',
  'Pm': 'Pm_3',
  'Pr': 'Pr_3',
  'Pt': 'Pt',
  'Pu': 'Pu',
  'Rb': 'Rb_sv',
  'Re': 'Re_pv',
  'Rh': 'Rh_pv',
  'Ru': 'Ru_pv',
  'S': 'S',
  'Sb': 'Sb',
  'Sc': 'Sc_sv',
  'Se': 'Se',
  'Si': 'Si',
  'Sm': 'Sm_3',
  'Sn': 'Sn_d',
  'Sr': 'Sr_sv',
  'Ta': 'Ta_pv',
  'Tb': 'Tb_3',
  'Tc': 'Tc_pv',
  'Te': 'Te',
  'Th': 'Th',
  'Ti': 'Ti_pv',
  'Tl': 'Tl_d',
  'Tm': 'Tm_3',
  'U': 'U',
  'V': 'V_pv',
  'W': 'W_pv',
  'Xe': 'Xe',
  'Y': 'Y_sv',
  'Yb': 'Yb_2',
  'Zn': 'Zn',
  'Zr': 'Zr_sv'}}

filename = input('Choose file name for vasp default parameters(default: vasp.yml): ')
filename = filename if filename else 'vasp.yml' 

with open(os.path.join(filepath,filename),'w+') as f:
    yaml.dump(vasp_default,f) 

print(f'Vasp default file {filename} written in {filepath}')  













