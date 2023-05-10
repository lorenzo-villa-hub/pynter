#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:49:47 2023

@author: villa
"""
import os
import os.path as op
from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.outputs import Vasprun

from pynter.vasp.jobs import VaspJob
from pynter.vasp.__init__ import load_vasp_default

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/vasp/tests/test_files')

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

def test_vaspjob_from_directory():
    path = op.join(test_files_path,'Si-BS')
    j = VaspJob.from_directory(path,load_outputs=True)
    assert j.inputs.as_dict() == VaspInput.from_directory(path).as_dict()
    assert j.outputs['Vasprun'].as_dict() == Vasprun(op.join(path,'vasprun.xml')).as_dict()
    assert j.job_settings == job_settings
    assert j.job_settings == job_settings    
    assert j.final_energy == -11.00288193
    assert j.charge == 0
    
def test_vaspjob_from_json():
    j = VaspJob.from_json(op.join(test_files_path,'Si-BS.json'))
    assert j.job_settings == job_settings    
    assert j.final_energy == -11.00288193
    assert j.charge == 0
    keys = load_vasp_default()['computed_entry_default']
    assert list(j.computed_entry.data.keys()) == keys
    
    