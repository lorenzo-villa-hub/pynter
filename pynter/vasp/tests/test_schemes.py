#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 12:04:33 2023

@author: villa
"""

import os.path as op

from pynter.vasp.jobs import VaspJob
from pynter.vasp.schemes import InputSets, Schemes, AdvancedSchemes
from pynter.data.datasets import Dataset

from pynter.testing.core import PynterTest
from pynter.testing.data import DatasetTest
from pynter.testing.vasp import VaspJobTest


class TestVaspSchemes(PynterTest):

    def test_get_vaspjob(self):
        j = InputSets(self.test_files_path,structure=self.structure,incar_settings=self.incar_settings,
                      name='Si_input_sets',job_settings=self.job_settings).get_vaspjob(setname='job')
        j_test = VaspJob.from_json(op.join(self.test_files_path,'vaspjob_Si_input_sets.json'))
        
        VaspJobTest().assert_inputs_equal(j,j_test,include_incar=False)
        return
        
    def test_scheme(self):
        ds_test = Dataset.from_directory(op.join(self.test_files_path,'Si_HSE_rel_gamma'))
        schemes = Schemes(self.test_files_path,structure=self.structure,incar_settings=self.incar_settings,
                      name='Si_schemes',job_settings=self.job_settings)
        jobs = schemes.hse_relaxation_gamma_extended()
        ds = Dataset(jobs)
        DatasetTest().assert_jobs_equal(ds.jobs, ds_test.jobs)
        return
        
    def test_advanced_scheme(self):
        ds_test = Dataset.from_directory(op.join(self.test_files_path,'Vac_Si_adv_schemes_inputs'))
        schemes = AdvancedSchemes(self.test_files_path,structure=self.structure,incar_settings=self.incar_settings,
                      name='Si_adv_schemes',job_settings=self.job_settings)
        jobs = schemes.vacancies_pbe_relaxation({'Si':[-1,0,1]},supercell_size=3)
        ds = Dataset(jobs)
        
        DatasetTest().assert_dataset_equal(ds, ds_test)
        return
        

    