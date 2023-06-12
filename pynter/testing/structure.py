#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""

import unittest
from numpy.testing import assert_array_almost_equal


class StructureTest(unittest.TestCase):
    
    def assertSiteEqual(self,site1,site2,decimal=3):
        assert_array_almost_equal(site1.frac_coords,site2.frac_coords,decimal=decimal)
        
    def assertLatticeEqual(self,lattice1,lattice2):
        assert_array_almost_equal(lattice1,lattice2)
        
    def assertStructureEqual(self,structure1,structure2):
        if len(structure1) != len(structure2):
            raise AssertionError
        else:
            for i in len(structure1):
                self.assertSiteEqual(structure1[i], structure2[i])
        
            

