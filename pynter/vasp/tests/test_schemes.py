#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 12:04:33 2023

@author: villa
"""

import os.path as op

from pymatgen.io.vasp.inputs import Kpoints

from pynter.vasp.jobs import VaspJob
from pynter.vasp.schemes import DefaultInputs, InputSets, Schemes, AdvancedSchemes
from pynter.data.datasets import Dataset

from pynter.testing.core import PynterTest
from pynter.testing.data import DatasetTest
from pynter.testing.vasp import VaspJobTest
from pynter.testing.vasp import VaspInputsTest


class TestDefaultInputs(PynterTest):

    def test_get_incar_default(self):
        di = DefaultInputs(structure=self.structure)
        incar_dict = di.get_incar_default()
        assert incar_dict["SYSTEM"] == "Si"
        assert incar_dict["ISYM"] == 2
    
    def test_get_kpoints_default(self):    
        di = DefaultInputs(structure=self.structure)
        kpoints = di.get_kpoints_default(kppa=1000)
        assert kpoints.kpts == [(8, 8, 8)]
        assert kpoints.kpts_shift == (0, 0, 0)
        
    def test_get_kpoints_bs_default(self):
        di = DefaultInputs(structure=self.structure)
        kpoints = di.get_kpoints_bs_default(divisions=10,hybrid_mode=False,kppa=1000)
        kpoints = Kpoints.from_str(str(kpoints)) # avoid pymatgen inconsistencies
        kpoints_file = Kpoints.from_file(op.join(self.test_files_path,'KPOINTS_bs_Si'))
        VaspInputsTest().assert_Kpoints_equal(kpoints, kpoints_file)
        
        kpoints = di.get_kpoints_bs_default(divisions=10,hybrid_mode=True,kppa=1000)
        kpoints = Kpoints.from_str(str(kpoints)) # avoid pymatgen inconsistencies
        kpoints_file = Kpoints.from_file(op.join(self.test_files_path,'KPOINTS_bs_hybrid_Si'))
        VaspInputsTest().assert_Kpoints_equal(kpoints, kpoints_file)
    
    
    def test_get_potcar(self):
        di = DefaultInputs(structure=self.structure)
        potcar = di.get_potcar(potcar_functional='PBE')
        assert len(potcar) == 1
        assert potcar[0].symbol == "Si"

class TestVaspSchemes(PynterTest):

    def test_get_vaspjob(self):
        j = InputSets(self.test_files_path,structure=self.structure,incar_settings=self.incar_settings,
                      name='Si_input_sets',job_settings=self.job_settings).get_vaspjob(setname='job')
        j_test = VaspJob.from_directory(op.join(self.test_files_path,'Si-input-sets'))
        
        VaspJobTest().assert_inputs_equal(j,j_test,include_incar=False)
        return
        
    def test_scheme(self):
        ds_test = Dataset.from_directory(op.join(self.test_files_path,'Si_HSE_rel_gamma'))
        job_settings = self.job_settings
        job_settings['slurm']['time'] = '24:00:00'
        schemes = Schemes(self.test_files_path,structure=self.structure,incar_settings=self.incar_settings,
                      name='Si_schemes',job_settings=job_settings)
        jobs = schemes.hse_relaxation_gamma_extended()       
        ds = Dataset(jobs)
        ds[3].job_settings['add_stop_array'] = True # stop array in job.sh file
        ds[7].job_settings['add_stop_array'] = True 
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
        

    