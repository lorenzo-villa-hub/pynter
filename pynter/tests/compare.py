#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""

import unittest
from numpy.testing import assert_array_almost_equal


class CompareStructures(unittest.TestCase):
    
    def compare_sites(self,site1,site2,decimal=3):
        assert_array_almost_equal(site1.frac_coords,site2.frac_coords,decimal=decimal)
        
    def compare_lattice(self,lattice1,lattice2):
        self.assertEqual(lattice1,lattice2)
        
    def compare(self,structure1,structure2):
        self.assertEqual(structure1.as_dict(),structure2.as_dict())
        
            

