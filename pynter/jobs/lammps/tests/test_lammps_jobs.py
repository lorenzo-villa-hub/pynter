#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 18:46:19 2025

@author: lorenzo
"""
import os.path as op
from pymatgen.core.structure import Structure

from pynter.jobs.lammps.lammps_jobs import LammpsJob

from pynter.testing.core import PynterTest

class TestLammpsJob(PynterTest):
            
    def test_from_directory_IN(self):
        job = LammpsJob.from_directory(op.join(self.test_files_path,'Si-IN'))
        assert job.job_settings['script_lines'][-1] == 'srun lmp -in input.in'
        assert job.inp.get_args('pair_style') == 'sw'
        assert job.formula == 'Si2'
        assert job.is_converged == None
        
    def test_from_directory_CG(self):
        job = LammpsJob.from_directory(op.join(self.test_files_path,'NBT-AL-more-fixes-CG'))
        assert type(job.initial_structure) == Structure
        assert job.formula == 'Na4 Ti8 Bi4 O24'
        assert job.inp.get_args('pair_style') == 'pace/extrapolation'
        assert all(key in job.outputs.keys() for key in ['log','convergence'])
        assert job.is_converged == True
        assert job.log[2]['TotEng'][3] == -306.17886
                
    def test_from_directory_FAILED(self):
        job = LammpsJob.from_directory(op.join(self.test_files_path,'STO-mod-FAILED'),
                                       job_script_filename='modified.sh',
                                       atom_style='charge',
                                       input_filename='modified.in',
                                       data_filename='modified.data',
                                       log_filename='modified.log.lammps')
        assert job.input_filename == 'modified.in'
        assert job.data_filename == 'modified.data'
        assert job.log_filename == 'modified.log.lammps'
        assert job.formula == 'Sr26 Ti27 O81'
        assert 'convergence' in job.outputs.keys()
        assert job.is_converged == False
        assert job.log == False
