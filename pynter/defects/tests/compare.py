#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""

import unittest

from pynter.tests.compare import CompareStructures

class CompareDefects(unittest.TestCase):
    
    def compare(self,defect1,defect2):
        CompareStructures().compare_sites(defect1.site,defect2.site)
        CompareStructures().compare(defect1.bulk_structure,defect2.bulk_structure)
        self.assertEqual(defect1.charge,defect2.charge)
        self.assertEqual(defect1.multiplicity,defect2.multiplicity)
        self.assertEqual(defect1.label,defect2.label)
        return
        
class CompareEntries(unittest.TestCase):
    
    def compare(self,entry1,entry2):
        CompareDefects().compare(entry1.defect,entry2.defect)
        CompareStructures().compare(entry1.bulk_structure,entry2.bulk_structure)
        self.assertEqual(entry1.energy_diff,entry2.energy_diff)
        self.assertEqual(entry1.corrections,entry2.corrections)
        self.assertEqual(entry1.data,entry2.data)
        self.assertEqual(entry1.label,entry2.label)
            
        
            

