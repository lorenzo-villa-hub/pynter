#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""

import unittest
from numpy.testing import assert_allclose

from pymatgen.io.vasp.inputs import Kpoints, Poscar

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
        kpoints1 = Kpoints.from_string(str(kpoints1))
        kpoints2 = Kpoints.from_string(str(kpoints2))
        assert_allclose(kpoints1.kpts, kpoints2.kpts)
        self.assertEqual(kpoints1.style,kpoints2.style)
        self.assertEqual(kpoints1.kpts_shift,kpoints2.kpts_shift)
        self.assertEqual(kpoints1.kpts_weights,kpoints2.kpts_weights)
        
    def compare_poscar(self,poscar1,poscar2):
        """
        Compare Poscar objects
        """
        poscar1 = Poscar.from_string(str(poscar1))
        poscar2 = Poscar.from_string(str(poscar2))
        self.assertEqual(poscar1.get_string(significant_figures=2),poscar2.get_string(significant_figures=2))
        
    def compare_potcar(self,potcar1,potcar2):
        """
        Compare Potcar objects
        """
        if len(potcar1) != len(potcar2):
            raise AssertionError
        else:
            for i in range(0,len(potcar1)):
                self.assertEqual(potcar1[i].element, potcar2[i].element)
                self.assertEqual(potcar1[i].functional, potcar2[i].functional)
     
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
            
            
            
            

