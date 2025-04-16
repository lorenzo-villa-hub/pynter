#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 17:52:51 2023

@author: lorenzo
"""
from numpy.testing import assert_allclose
from pymatgen.io.vasp.inputs import Kpoints, Poscar, Incar
from pynter.testing.core import PynterTest


class VaspInputsTest(PynterTest):
    """
    Provides methods to test pymatgen.io.vasp.inputs objects
    """
    def assert_Incar_equal(self, incar1, incar2, system_only=False):
        """
        Compare Incar objects. Use system_only when default parameters are not known.
        """
        for inc in [incar1,incar2]:
            if type(inc) != Incar:
                inc = Incar(inc)
            for k,v in inc.items():
                inc.__setitem__(k,v)
        if system_only:
            self.assertEqual(incar1['SYSTEM'], incar2['SYSTEM'])
        else:
            assert incar1.diff(incar2)['Different'] == {}

    def assert_Kpoints_equal(self,kpoints1,kpoints2):
        """
        Compare Kpoint objects, compares the strings to avoid inconsistencies
        which can happen when importing/exporting to dict
        """
        kpoints1 = Kpoints.from_str(str(kpoints1))
        kpoints2 = Kpoints.from_str(str(kpoints2))
        assert_allclose(kpoints1.kpts, kpoints2.kpts,atol=1e-03)
        self.assertEqual(kpoints1.style,kpoints2.style)
        self.assertEqual(kpoints1.kpts_shift,kpoints2.kpts_shift)
        self.assertEqual(kpoints1.kpts_weights,kpoints2.kpts_weights)
        
    def assert_Poscar_equal(self,poscar1,poscar2):
        """
        Compare Poscar objects
        """
        poscar1 = Poscar.from_str(str(poscar1))
        poscar2 = Poscar.from_str(str(poscar2))
        self.assertEqual(poscar1.get_string(significant_figures=4),poscar2.get_string(significant_figures=4))
        
    def assert_Potcar_equal(self,potcar1,potcar2):
        """
        Compare Potcar objects
        """
        if len(potcar1) != len(potcar2):
            raise AssertionError
        else:
            for i in range(0,len(potcar1)):
                self.assertEqual(potcar1[i].element, potcar2[i].element)
                self.assertEqual(potcar1[i].functional, potcar2[i].functional)
     
    def assert_VaspInput_equal(self,input1,input2,include_incar=True,system_only=False):
        if include_incar:
            self.assert_Incar_equal(input1['INCAR'],input2['INCAR'],system_only)
        self.assert_Kpoints_equal(input1['KPOINTS'],input2['KPOINTS'])
        self.assert_Poscar_equal(input1['POSCAR'],input2['POSCAR'])
        self.assert_Potcar_equal(input1['POTCAR'],input2['POTCAR'])
        
        
class VaspOutputsTest(PynterTest):
    """
    Provides methods to test pymatgen.io.vasp.outputs objects and pymatgen ComputedEntry
    """
    def assert_ComputedEntry_equal(self,entry1,entry2):
        attributes = ['correction','correction_uncertainty','energy_adjustments',
                      'energy','uncorrected_energy']
        for attr in attributes:
            value1, value2 = getattr(entry1,attr), getattr(entry2,attr)
            self.assert_object_almost_equal(value1,value2)   
        self.assertEqual(entry1.data.keys(),entry2.data.keys())
            
    
    def assert_Vasprun_equal(self,vasprun1,vasprun2):
        self.assert_object_almost_equal(vasprun1.as_dict(),vasprun2.as_dict(),rtol=0.001)
        

        
        
        
        
        
        