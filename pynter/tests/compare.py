#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""

import unittest
#from pymatgen.io.vasp.inputs import Incar,Kpoints,Poscar,Potcar, VaspInput
from pynter.vasp.jobs import VaspJob


class CompareVaspInputs(unittest.TestCase):
    
    def compare_incar(self, incar1, incar2, system_only=False):
        """
        Compare Incar objects. Use system_only when default parameters are not known.
        """
        if system_only:
            self.assertEqual(incar1['SYSTEM'], incar2['SYSTEM'])
        else:
            self.assertEqual(incar1,incar2)

    def compare_kpoints(self,kpoints1,kpoints2):
        """
        Compare Kpoint objects, compares the strings to avoid inconsistencies
        which can happen when importing/exporting to dict
        """
        self.assertEqual(str(kpoints1), str(kpoints2))
        
    def compare_poscar(self,poscar1,poscar2):
        """
        Compare Poscar objects
        """
        self.assertEqual(poscar1.as_dict(),poscar2.as_dict())
        
    def compare_potcar(self,potcar1,potcar2):
        """
        Compare Potcar objects
        """
        self.assertEqual(potcar1.as_dict(),potcar2.as_dict())
        
    def compare_vaspinput(self,input1,input2,include_incar=True,system_only=False):
        if include_incar:
            self.compare_incar(input1['INCAR'],input2['INCAR'],system_only)
        self.compare_kpoints(input1['KPOINTS'],input2['KPOINTS'])
        self.compare_poscar(input1['POSCAR'],input2['POSCAR'])
        self.compare_potcar(input1['POTCAR'],input2['POTCAR'])
        
        
class CompareVaspOutputs(unittest.TestCase):
    
    def compare_computed_entry(self,entry1,entry2):
        self.assertEqual(entry1.as_dict(),entry2.as_dict())
    
    def compare_vasprun(self,vasprun1,vasprun2):
        self.assertEqual(vasprun1.as_dict(),vasprun2.as_dict())


class CompareJobs(unittest.TestCase):
    
    def __init__(self,jobclass='VaspJob'):
        self.jobclass = jobclass
        if self.jobclass == 'VaspJob':
            self.cls = CompareVaspJobs
    
    def compare_inputs(self,job1,job2,**kwargs):
        self.cls().compare_inputs(**kwargs)
        
    def compare_outputs(self,job1,job2,**kwargs):
        self.cls().compare_outputs(**kwargs)
        
    def compare_settings(self,job1,job2):
        self.assertEqual(job1.job_settings,job2.job_settings)
        self.assertEqual(job1.name,job2.name)
        
    def compare(self,job1,job2,kwargs_input={},kwargs_output={}):
        self.compare_inputs(job1,job2,**kwargs_input)
        if job1.outputs or job2.outputs:
            self.compare_outputs(job1,job2)
        self.compare_settings(job1,job2)
        

class CompareVaspJobs(CompareJobs):
    
    def compare_inputs(self,job1,job2,include_incar=True,system_only=False):
        CompareVaspInputs().compare_vaspinput(job1.inputs,job2.inputs,include_incar,system_only)

    def compare_outputs(self,job1,job2):
        CompareVaspOutputs().compare_computed_entry(job1.computed_entry,job2.computed_entry)
        if "Vasprun" in [job1.outputs.keys(),job2.outputs.keys()]:
            CompareVaspOutputs().compare_vasprun(job1.vasprun, job2.vasprun)
            
            
class CompareDatasets(unittest.TestCase):
    
    def compare_jobs(self,ds1,ds2,kwargs_input={},kwargs_output={}):
        ds1.sort_jobs()
        ds2.sort_jobs()
        for i in range(len(ds1.jobs)):
            job1,job2 = ds1.jobs[i], ds2.jobs[i]
            CompareJobs().compare(job1, job2,kwargs_input,kwargs_output)
            
            

