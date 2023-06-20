#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:39:37 2023

@author: villa
"""
import yaml
from setuptools import setup, find_namespace_packages
import atexit
import os

setup(
    name='pynter-defects',
    version='1.0.3',
    author='Lorenzo Villa',
    description='Tools for atomistic calculations, provides features for point-defect calculations with VASP',
    packages=find_namespace_packages(exclude=["pynter.*.tests", "pynter.*.*.tests"]),
    install_requires=[
        'ase',
        'pymatgen>=2023.3.23',
        'pymatgen-analysis-defects>=2023.3.25',
        'pymatgen-db',
        'PyYAML',
        'schedule'        
    ],
    extra_requires={'test':'pytest'},
    entry_points={
    "console_scripts": [
        "pynter = pynter.cli.main:main"]}
)


homedir = os.getenv("HOME")
settings_path = os.path.join(homedir,'.pynter')
config_filename = 'config.yml'
vasp_config_filename = 'vasp.yml'

def post_install():
    
    if not os.path.exists(os.path.join(settings_path,config_filename)):
        set_default_config()
    if not os.path.exists(os.path.join(settings_path,vasp_config_filename)):
        set_default_vasp_config()
    return


def set_default_config():
    config = {        
        'HPC': 
              {'hostname': None,
              'localdir': os.path.join(homedir,'pynter-data'),
              'workdir': None},
        'API_KEY': None,
        'job_settings': 
              {'project_id': '',
              'name': 'no_name',
              'array_size': None,
              'email': None,
              'nodes': 1,
              'cores_per_node': 1,
              'output_filename': 'out.%j',
              'error_filename': 'err.%j',
              'timelimit': '01:00:00',
              'memory_per_cpu': 3500,
              'partition': None,
              'processor': None,
              'modules': None,
              'path_exe': '',
              'add_stop_array': True,
              'add_automation': None,
              'add_lines_header': None,
              'add_lines_body': None,
              'filename': 'job.sh'},
        'dbconfig': 
              {'vasp': 
                   {'collection': 'vasp',
                   'database': 'pynter',
                   'host': '127.0.0.1',
                   'password': None,
                   'port': 27017,
                   'user': None}}}
    
    with open(os.path.join(settings_path,config_filename),'w+') as f:
        yaml.dump(config,f) 
    return


def set_default_vasp_config():
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
    
    with open(os.path.join(settings_path,vasp_config_filename),'w+') as f:
        yaml.dump(vasp_default,f) 
    return

atexit.register(post_install)