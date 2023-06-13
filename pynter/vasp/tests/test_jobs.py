#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:49:47 2023

@author: villa
"""
import os.path as op
import json

from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.outputs import Vasprun

from pynter.vasp.jobs import VaspJob
from pynter.vasp.__init__ import load_vasp_default
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

        assert j.job_settings['name'] == 'Si-BS_PBE-el-str_3'  
        self.assert_all_close(j.final_energy,-11.00288193)
        self.assert_all_close(j.charge,0)
    
    def test_vaspjob_from_json(self):
        j = VaspJob.from_json(op.join(self.test_files_path,'Si-BS.json'))
        assert j.job_settings == self.job_settings    
        assert j.final_energy == -11.00288193
        assert j.charge == 0
        keys = load_vasp_default()['computed_entry_default']
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
        

    