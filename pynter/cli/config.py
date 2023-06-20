#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:24:40 2023

@author: villa
"""

import os
import yaml


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
        email = input('email address (Default None): ')
        path_exe = input('Path of default executable: ')
        job_script_filename = input('Job script filename (default job.sh): ')
        
        API_KEY = input('\nMaterials Project API_KEY (default None):')
        
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
            'API_KEY': API_KEY,
            'job_settings': 
                  {'project_id': project_id,
                  'name': 'no_name',
                  'array_size': None,
                  'email': email,
                  'nodes': 4,
                  'cores_per_node': 96,
                  'output_filename': 'out.%j',
                  'error_filename': 'err.%j',
                  'timelimit': '24:00:00',
                  'memory_per_cpu': 3500,
                  'partition': None,
                  'processor': None,
                  'modules': ['intel/2020.2','intelmpi/2020.2','fftw/3.3.8'],
                  'path_exe': path_exe,
                  'add_stop_array': True,
                  'add_automation': None,
                  'add_lines_header': ['I_MPI_PMI_LIBRARY=/opt/slurm/current/lib/libpmi2.so'],
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
     
     'computed_entry_default': ['final_energy','structures','eigenvalue_band_properties',
                                'parameters','actual_kpoints','ionic_steps'],
     
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
      'Na': 'Na_pv',
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
    
    # set vasp.yml
    if exclude != filename_vasp:
    
        filename = filename_vasp
        
        with open(os.path.join(filepath,filename),'w+') as f:
            yaml.dump(vasp_default,f) 
        
        print(f'Vasp default file {filename} written in {filepath}')  













