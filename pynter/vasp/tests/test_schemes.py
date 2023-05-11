#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 12:04:33 2023

@author: villa
"""

import os
import os.path as op
import json
from monty.json import MontyEncoder
from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.outputs import Vasprun

from pynter.tests.__init__ import get_structure_Si, get_job_settings, get_incar_settings
from pynter.vasp.jobs import VaspJob
from pynter.vasp.schemes import InputSets, Schemes
from pynter.data.datasets import Dataset
from pynter.tests.compare import CompareVaspJobs

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/vasp/tests/test_files')

structure = get_structure_Si()
job_settings = get_job_settings()
incar_settings = get_incar_settings()

def test_get_vaspjob():
    j = InputSets(test_files_path,structure=structure,incar_settings=incar_settings,
                  name='Si_input_sets',job_settings=job_settings).get_vaspjob(setname='job')
    j_test = VaspJob.from_json(op.join(test_files_path,'Si_input_sets.json'))
    CompareVaspJobs().compare_inputs(j,j_test,include_incar=False)
    
    
def test_scheme():
    ds_test = Dataset.from_json(op.join(test_files_path,'Si_HSE_rel_gamma_inputs.json'))
    schemes = Schemes(test_files_path,structure=structure,incar_settings=incar_settings,
                  name='Si_input_sets',job_settings=job_settings)
    jobs = schemes.hse_ionic_rel_gamma()
    ds = Dataset(jobs)
    
    