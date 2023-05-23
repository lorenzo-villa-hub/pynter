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
from pynter.vasp.schemes import InputSets, Schemes, AdvancedSchemes
from pynter.data.datasets import Dataset
from pynter.vasp.tests.compare import CompareVaspJobs
from pynter.data.tests.compare import CompareDatasets

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/vasp/tests/test_files')

structure = get_structure_Si()
job_settings = get_job_settings()
incar_settings = get_incar_settings()


def test_get_vaspjob():
    j = InputSets(test_files_path,structure=structure,incar_settings=incar_settings,
                  name='Si_input_sets',job_settings=job_settings).get_vaspjob(setname='job')
    j_test = VaspJob.from_json(op.join(test_files_path,'vaspjob_Si_input_sets.json'))
    CompareVaspJobs().compare_inputs(j,j_test,include_incar=False)
    
    
def test_scheme():
    ds_test = Dataset.from_json(op.join(test_files_path,'ds_Si_HSE_rel_gamma_inputs.json'))
    schemes = Schemes(test_files_path,structure=structure,incar_settings=incar_settings,
                  name='Si_schemes',job_settings=job_settings)
    jobs = schemes.hse_relaxation_gamma_extended()
    ds = Dataset(jobs)
    CompareDatasets().compare_jobs(ds, ds_test)
    
    
def test_advanced_scheme():
    ds_test = Dataset.from_json(op.join(test_files_path,'ds_Si_vacancies_pbe_relaxation.json'))
    schemes = AdvancedSchemes(test_files_path,structure=structure,incar_settings=incar_settings,
                  name='Si_adv_schemes',job_settings=job_settings)
    jobs = schemes.vacancies_pbe_relaxation({'Si':[-1,0,1]},supercell_size=3)
    ds = Dataset(jobs)
    CompareDatasets().compare_jobs(ds, ds_test)
        

    