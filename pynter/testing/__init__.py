#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:13:48 2023

@author: villa
"""

from pymatgen.core.structure import Structure

def get_structure_Si():
    structure = Structure.from_dict(
    {'@module': 'pymatgen.core.structure',
     '@class': 'Structure',
     'charge': 0,
     'lattice': {'matrix': [[3.32548851, 0.0, 1.91997169],
       [1.10849617, 3.13530064, 1.91997169],
       [0.0, 0.0, 3.83994338]],
      'pbc': (True, True, True),
      'a': 3.839943374653261,
      'b': 3.839943378813096,
      'c': 3.83994338,
      'alpha': 59.99999998977525,
      'beta': 59.99999995393976,
      'gamma': 60.00000000512866,
      'volume': 40.036809671145996},
     'sites': [{'species': [{'element': 'Si', 'occu': 1}],
       'abc': [0.875, 0.875, 0.875],
       'xyz': [3.879736595, 2.74338806, 6.719900914999999],
       'label': 'Si',
       'properties': {}},
      {'species': [{'element': 'Si', 'occu': 1}],
       'abc': [0.125, 0.125, 0.125],
       'xyz': [0.554248085, 0.39191258, 0.959985845],
       'label': 'Si',
       'properties': {}}]}
        )
    return structure



def get_incar_settings():
    incar_settings = {'ALGO': 'Normal',
     'AMIX': 0.2,
     'EDIFF': 1e-06,
     'EDIFFG': -0.05,
     'ENCUT': 550,
     'IBRION': 2,
     'ISIF': 2,
     'ISMEAR': 0,
     'ISPIN': 1,
     'KPAR': 4,
     'LCHARG': '.TRUE.',
     'LORBIT': 10,
     'LREAL': '.FALSE.',
     'LWAVE': '.TRUE.',
     'NELM': 200,
     'NSW': 0,
     'SIGMA': 0.05,
     'SYSTEM': 'Si',
     '#### Default PBE: system': 'Si',
     'ISYM': 2}
    return incar_settings



def get_job_settings():
    job_settings = {
     'add_automation': 'automation_vasp.py --contcar --chgcar --wavecar --check-kpoints --error-check',
     'add_lines_body': None,
     'add_lines_header': None,
     'add_stop_array': True,
     'array_size': None,
     'cores_per_node': 24,
     'email': 'test@pynter.com',
     'error_filename': 'err.%j',
     'filename': 'job.sh',
     'memory_per_cpu': 2400,
     'modules': ['intel/2019.2', 'intel/2019.3', 'intelmpi/2019.3', 'fftw/3.3.8'],
     'name': 'Si-BS_PBE-el-str_3',
     'nodes': 1,
     'output_filename': 'out.%j',
     'partition': 'deflt',
     'path_exe': '/home/vasp-5-3-3',
     'processor': 'avx2',
     'project_id': 'project0000',
     'timelimit': '00:30:00'}
    return job_settings