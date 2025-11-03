#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 12:04:33 2023

@author: villa
"""

import os.path as op

from pymatgen.io.vasp.inputs import Kpoints

from pynter.hpc.slurm import JobSettings

from pynter.jobs.vasp.vasp_jobs import VaspJob
from pynter.jobs.datasets import Dataset

from pynter.vasp.schemes.core import DefaultInputs, DefaultJobSettings, InputSets
from pynter.vasp.schemes.relaxation import RelaxationSchemes
from pynter.vasp.schemes.defects import DefectSchemes


from pynter.testing.core import PynterTest
from pynter.testing.jobs import DatasetTest, VaspJobTest
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


class TestDefaultJobSettings(PynterTest):
    
    def setUp(self):
        self.default = DefaultJobSettings(
            vasp_exe = '/test/vasp',
            vasp_sbatch={'ntasks':384},
            modules=['module1','module2'],
            lines_before_srun=['line1 before','line2 before'],
            lines_after_srun=['line1 after','line2 after']
            )
    
    def test_vasp_script_lines(self):        
        updated_js = self.default.get_updated_job_settings(JobSettings())
        lines = updated_js['script_lines']

        self.assertIn('module load module1',lines)
        self.assertIn('module load module2',lines)
        assert lines.index('module load module1') < lines.index('module load module2')
        self.assertIn('line1 before',lines)
        self.assertIn('line2 before',lines)
        assert lines.index('line1 before') < lines.index('line2 before')
        self.assertIn('line1 after',lines)
        self.assertIn('line2 after',lines)
        assert lines.index('line1 after') < lines.index('line2 after')
        
    def test_duplicates_flexibility(self):
        original_js = self.default.get_updated_job_settings(JobSettings())
        updated_js = DefaultJobSettings(vasp_exe='/test/vasp').get_updated_job_settings(original_js)
        def assert_no_duplicates(target_string,lines):
            count = 0
            for line in lines:
                if target_string in line:
                    count += 1
            assert count <= 1
        strings_to_check = ['module load module1','line1 before','srun','line1 after']
        for string in strings_to_check:
            assert_no_duplicates(string, updated_js['script_lines'])
        


class TestVaspSchemes(PynterTest):

    def test_get_vaspjob(self):
        j = InputSets(self.test_files_path,structure=self.structure,incar_settings=self.incar_settings,
                      name='Si_input_sets',job_settings=self.job_settings).get_vaspjob(add_to_job_name='job')
        j_test = VaspJob.from_directory(op.join(self.test_files_path,'Si-input-sets'))
        
        VaspJobTest().assert_inputs_equal(j,j_test,include_incar=False)
        return
        
    def test_scheme(self):
        ds_test = Dataset.from_directory(op.join(self.test_files_path,'Si_HSE_rel_gamma'))
        job_settings = self.job_settings
        job_settings['sbatch']['time'] = '24:00:00'
        schemes = RelaxationSchemes(self.test_files_path,structure=self.structure,incar_settings=self.incar_settings,
                      name='Si_schemes',job_settings=job_settings)
        jobs = schemes.hse_relaxation_gamma_extended()       
        ds = Dataset(jobs)
        DatasetTest().assert_jobs_equal(ds.jobs, ds_test.jobs)
        return
        
    def test_defect_scheme(self):
        ds_test = Dataset.from_directory(op.join(self.test_files_path,'Vac_Si_adv_schemes_inputs'))
        schemes = DefectSchemes(self.test_files_path,structure=self.structure,incar_settings=self.incar_settings,
                      name='Si_adv_schemes',job_settings=self.job_settings)
        jobs = schemes.vacancies_pbe_relaxation({'Si':[-1,0,1]},supercell_size=3)
        ds = Dataset(jobs)
        
        DatasetTest().assert_dataset_equal(ds, ds_test)
        return
        

    