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
        
            

