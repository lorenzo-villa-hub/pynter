#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:49:47 2023

@author: villa
"""
import unittest
import os
import os.path as op
import json
from monty.json import MontyEncoder
from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.outputs import Vasprun

from pynter.vasp.jobs import VaspJob
from pynter.vasp.__init__ import load_vasp_default
from pynter.tests.__init__ import get_job_settings
from pynter.vasp.tests.compare import CompareVaspInputs, CompareVaspJobs

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/vasp/tests/test_files')

job_settings = get_job_settings()

    
def test_vaspjob_from_directory():
    path = op.join(test_files_path,'Si-BS')
    j = VaspJob.from_directory(path,load_outputs=True)
    assert j.inputs.as_dict() == VaspInput.from_directory(path).as_dict()
    assert j.outputs['Vasprun'].as_dict() == Vasprun(op.join(path,'vasprun.xml')).as_dict()
    assert j.job_settings['name'] == job_settings['name']  
    assert j.final_energy == -11.00288193
    assert j.charge == 0

def test_vaspjob_from_json():
    j = VaspJob.from_json(op.join(test_files_path,'Si-BS.json'))
    assert j.job_settings == job_settings    
    assert j.final_energy == -11.00288193
    assert j.charge == 0
    keys = load_vasp_default()['computed_entry_default']
    assert list(j.computed_entry.data.keys()) == keys
    
def test_vaspjob_with_band_structure():
    path = op.join(test_files_path,'Si-BS')
    j = VaspJob.from_directory(path,load_outputs=True)
    j.get_output_properties(get_band_structure=True)
    assert j.band_structure.as_dict() == Vasprun(op.join(path,'vasprun.xml')).get_band_structure().as_dict()
    
def test_vaspjob_to_json_from_json():
    path = op.join(test_files_path,'Si-BS')
    j = VaspJob.from_directory(path,load_outputs=True)
    j.get_output_properties(get_band_structure=True)
    d = j.as_dict(get_band_structure=True) 
    #if runs import is successfull - comparing dictionary directly is impossible because of pymatgen's inconsistencies
    j_from_json = VaspJob.from_json(json.dumps(d))
    CompareVaspJobs().compare(j, j_from_json)
    assert j_from_json.band_structure.as_dict() == j.band_structure.as_dict()

    