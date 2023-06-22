#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""
from pynter.testing.core import PynterTest
from pynter.testing.structure import StructureTest


class DefectTest(PynterTest):
    """
    Provides methods to test Defect objects
    """
    def assert_Defect_equal(self,defect1,defect2):
        StructureTest().assert_Site_equal(defect1.site,defect2.site)
        StructureTest().assert_Structure_equal(defect1.bulk_structure,defect2.bulk_structure)
        self.assertEqual(defect1.charge,defect2.charge)
        self.assertEqual(defect1.multiplicity,defect2.multiplicity)
        self.assertEqual(defect1.label,defect2.label)
        
        
class DefectEntryTest(PynterTest):
    """
    Provides methods to test DefectEntry objects
    """
    def assert_DefectEntry_equal(self,entry1,entry2):
        DefectTest().assert_Defect_equal(entry1.defect,entry2.defect)
        StructureTest().assert_Structure_equal(entry1.bulk_structure,entry2.bulk_structure)
        self.assertEqual(entry1.energy_diff,entry2.energy_diff)
        self.assertEqual(entry1.corrections,entry2.corrections)
        self.assertEqual(entry1.data,entry2.data)
        self.assertEqual(entry1.label,entry2.label)
        

            
        
            

