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
        self.path = op.join(self.test_files_path,'Si-BS')
        self.job = VaspJob.from_directory(self.path,load_outputs=True)

    def test_vaspjob_from_directory(self):
        path, j = self.path, self.job        
        VaspInputsTest().assert_VaspInput_equal(j.inputs,VaspInput.from_directory(path))
        VaspOutputsTest().assert_Vasprun_equal(j.outputs['Vasprun'], Vasprun(op.join(path,'vasprun.xml')))

        assert j.job_settings['job-name'] == 'Si-BS_PBE-el-str_3'  
        self.assert_all_close(j.final_energy,-11.00288193)
        self.assert_all_close(j.charge,0)
    
    def test_vaspjob_from_dict(self):
        j = VaspJob.from_dict(VaspJob.from_directory(self.path,load_outputs=True).as_dict())
        keys = SETTINGS['vasp']['computed_entry_default']
        assert list(j.computed_entry.data.keys()) == keys
        
    def test_vaspjob_with_band_structure(self):
        j = self.job        
        j.get_output_properties(get_band_structure=True)
        assert j.band_structure.__module__ == 'pymatgen.electronic_structure.bandstructure'
        
    def test_vaspjob_to_json_from_json(self):
        j = self.job        
        j.get_output_properties(get_band_structure=True)
        d = j.as_dict(get_band_structure=True) 
        #if runs import is successfull - comparing dictionary directly is impossible because of pymatgen's inconsistencies
        j_from_json = VaspJob.from_json(json.dumps(d))
        assert j.name == 'Si-BS_PBE-el-str_3'
        

class TestVaspNEBJob(PynterTest):
    
    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore')
        self.path = op.join(self.test_files_path,'NN-VO-NEB')
        self.job = VaspNEBJob.from_directory(self.path)        
    
    def test_vaspnebjob_from_directory(self):
        neba = self.job.neb_analysis
        r = np.array([0.        , 0.68575757, 1.36377546, 2.0248057 , 2.77804612,
               3.5241437 ])
        self.assert_all_close( neba.r , r)
        
        energies = np.array([-1228.22426719, -1228.12890034, -1227.95020108, -1227.90419223,
               -1228.07030217, -1228.19596946])
        self.assert_all_close( neba.energies , energies)
        
        forces = np.array([ 0.      , -0.27551 , -0.22473 ,  0.089532,  0.337205,  0.      ])
        self.assert_all_close( neba.forces , forces)
        
    def test_vaspnebjob_from_json(self):
        #if runs import is successfull - comparing dictionary directly is impossible because of pymatgen's inconsistencies
        j_from_json = VaspNEBJob.from_json(json.dumps(self.job.as_dict()))
        return
        
        
        