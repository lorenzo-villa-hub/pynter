#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:49:47 2023

@author: villa
"""
import os.path as op
import json
import numpy as np

from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.outputs import Vasprun

from pynter.vasp.jobs import VaspJob, VaspNEBJob
from pynter import SETTINGS
from pynter.testing.core import PynterTest
from pynter.testing.vasp import VaspInputsTest, VaspOutputsTest


class TestVaspJob(PynterTest):
    
    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore')

    def test_vaspjob_from_directory(self):
        path = op.join(self.test_files_path,'Si-BS')
        j = VaspJob.from_directory(path,load_outputs=True)
        
        VaspInputsTest().assert_VaspInput_equal(j.inputs,VaspInput.from_directory(path))

        VaspOutputsTest().assert_Vasprun_equal(j.outputs['Vasprun'], Vasprun(op.join(path,'vasprun.xml')))

        assert j.job_settings['job-name'] == 'Si-BS_PBE-el-str_3'  
        self.assert_all_close(j.final_energy,-11.00288193)
        self.assert_all_close(j.charge,0)
    
    def test_vaspjob_from_json(self):
        j = VaspJob.from_json(op.join(self.test_files_path,'Si-BS.json'))
        assert j.job_settings == self.job_settings    
        assert j.final_energy == -11.00288193
        assert j.charge == 0
        keys = SETTINGS['vasp']['computed_entry_default']
        assert list(j.computed_entry.data.keys()) == keys
        
    def test_vaspjob_with_band_structure(self):
        path = op.join(self.test_files_path,'Si-BS')
        j = VaspJob.from_directory(path,load_outputs=True)
        j.get_output_properties(get_band_structure=True)
        assert j.band_structure.__module__ == 'pymatgen.electronic_structure.bandstructure'
        
    def test_vaspjob_to_json_from_json(self):
        path = op.join(self.test_files_path,'Si-BS')
        j = VaspJob.from_directory(path,load_outputs=True)
        j.get_output_properties(get_band_structure=True)
        d = j.as_dict(get_band_structure=True) 
        #if runs import is successfull - comparing dictionary directly is impossible because of pymatgen's inconsistencies
        j_from_json = VaspJob.from_json(json.dumps(d))
        assert j.name == 'Si-BS_PBE-el-str_3'
        

class TestVaspNEBJob(PynterTest):
    
    def test_vaspnebjob_from_json(self):
        j = VaspNEBJob.from_json(self.get_testfile_path('vaspnebjob_NN_VO.json'))
        neba = j.neb_analysis
        r = np.array([0.        , 0.64591638, 1.29416299, 1.94391391, 2.58088582,
               3.2171647 , 3.85151573])
        self.assert_all_close( neba.r , r)
        
        energies = np.array([-611.37449618, -611.20389509, -610.79989855, -610.64320248,
       -610.78631166, -611.14814619, -611.37438141])
        self.assert_all_close( neba.energies , energies)
        
        forces = np.array([ 0.      , -0.572757, -0.500791, -0.004058,  0.447456,  0.638393,
        0.      ])
        self.assert_all_close( neba.forces , forces)
        
        