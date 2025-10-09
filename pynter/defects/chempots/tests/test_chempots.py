#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:02:03 2023

@author: villa
"""
import numpy as np
from pymatgen.core.periodic_table import Element

from pynter.defects.chempots.core import Chempots, chempot_ideal_gas

from pynter.testing.core import PynterTest
from pynter.testing.chempots import ChempotsTest


class TestChempots(PynterTest):
    
    def setUp(self):
        self.mu_refs = {'Na': -1.31, 'Nb': -10.1, 'O': -4.95}
        self.mu = Chempots({'Na':-2.26,'Nb':-6.15,'O':-1.92})
        self.mu_abs = Chempots({'Na':-3.57,'Nb':-16.25,'O':-6.87})
        

    def test_as_dict_from_dict(self):
        self.assert_object_almost_equal( self.mu.as_dict() , Chempots.from_dict(self.mu.as_dict()).as_dict() )
        
    def test_pmg_elements(self):
        mu_elements = {Element('Na'):-2.26,Element('Nb'):-6.15,Element('O'):-1.92}
        ChempotsTest().assert_Chempots_equal( self.mu , Chempots.from_pmg_elements(mu_elements) )
        self.assert_object_almost_equal( self.mu.to_pmg_elements() , mu_elements )
        
    def test_absolute_referenced(self):
        ChempotsTest().assert_Chempots_equal( self.mu.get_absolute(self.mu_refs) , self.mu_abs, rtol=1e-03 )
        ChempotsTest().assert_Chempots_equal( self.mu_abs.get_referenced(self.mu_refs) , self.mu, rtol=1e-03 )
        
        
def test_chempot_ideal_gas():
    actual = chempot_ideal_gas(0,300,0.2)
    desired = -0.0208035922029101
    np.testing.assert_allclose(actual, desired, atol=1e-04)
    actual = chempot_ideal_gas(0,1300,0.2)
    desired = -0.09014889954594378
    np.testing.assert_allclose(actual, desired, atol=1e-04)
    
    

    
    
    
    
    
    
    
    
    
    