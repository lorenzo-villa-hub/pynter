#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""

import unittest

from pynter.data.tests.compare import CompareJobs


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
    
    def compare_computed_entry(self,entry1,entry2,flexible=True):
        if flexible:
            attributes = ['correction','correction_uncertainty','energy_adjustments',
                          'energy','uncorrected_energy']
            for attr in attributes:
                value1, value2 = getattr(entry1,attr), getattr(entry2,attr)
                self.assertEqual(value1,value2)   
            self.assertEqual(entry1.data.keys(),entry2.data.keys())
        else:
            self.assertEqual(entry1.as_dict(),entry2.as_dict())
    
    def compare_vasprun(self,vasprun1,vasprun2):
        self.assertEqual(vasprun1.as_dict(),vasprun2.as_dict())



        

class CompareVaspJobs(CompareJobs,unittest.TestCase):
    
    def compare_inputs(self,job1,job2,include_incar=True,system_only=False):
        CompareVaspInputs().compare_vaspinput(job1.inputs,job2.inputs,include_incar,system_only)

    def compare_outputs(self,job1,job2,flexible=True):
        CompareVaspOutputs().compare_computed_entry(job1.computed_entry,job2.computed_entry,flexible)
        if "Vasprun" in [job1.outputs.keys(),job2.outputs.keys()]:
            CompareVaspOutputs().compare_vasprun(job1.vasprun, job2.vasprun)
            
            
            
            

