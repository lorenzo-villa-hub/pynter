#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 10:08:17 2025

@author: lorenzo
"""
import os

from pymatgen.analysis.phase_diagram import PhaseDiagram

from pynter.tools.utils import get_object_from_json
from pynter.defects.chempots.core import Chempots
from pynter.defects.chempots.reservoirs import Reservoirs


from pynter.testing.chempots import ChempotsTest, ReservoirsTest
from .test_chempots import TestChempots

class TestReservoirs(TestChempots):
        
    def setUp(self):
        super().setUp()
        res_dict = {'X':self.mu}
        pd = get_object_from_json(PhaseDiagram,self.get_testfile_path('PD_Na-Nb-O.json'))
        self.res = Reservoirs(res_dict,pd,are_chempots_delta=True)
    
    def test_mu_refs(self):
        ChempotsTest().assert_Chempots_equal( self.res.mu_refs , Chempots(self.mu_refs), rtol=1e-02 )
        
    def test_as_dict_from_dict(self):
        ReservoirsTest().assert_Reservoirs_equal( self.res , Reservoirs.from_dict(self.res.as_dict()) ,rtol=1e-07)
        
    def test_to_json_from_json(self):
        try:
            self.res.to_json('reservoirs_temp.json')
            ReservoirsTest().assert_Reservoirs_equal( self.res , 
                                        Reservoirs.from_json('reservoirs_temp.json'), rtol=1e-02)
        finally:
            os.remove('reservoirs_temp.json')
    
    
    def test_absolute_referenced(self):
        self.res.set_to_absolute()
        ReservoirsTest().assert_Reservoirs_equal( self.res , Reservoirs({'X':self.mu_abs}) , check_reference=False, rtol=1e-03 )
        self.res.set_to_referenced()
        ReservoirsTest().assert_Reservoirs_equal( self.res , Reservoirs({'X':self.mu}), check_reference=False, rtol=1e-03)
        res_filtered = Reservoirs({'X':Chempots({'Na':-2.26,'O':-1.92})})
        ReservoirsTest().assert_Reservoirs_equal( 
            self.res.filter_reservoirs(elements=['Na','O']), res_filtered , check_reference=False, rtol=1e-03)