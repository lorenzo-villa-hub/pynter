#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 10:18:19 2025

@author: lorenzo
"""

from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram

from pynter.tools.utils import get_object_from_json
from pynter.defects.chempots.core import Chempots
from pynter.defects.chempots.reservoirs import Reservoirs
from pynter.defects.chempots.phase_diagram import PDHandler

from pynter.testing.core import PynterTest
from pynter.testing.chempots import  ReservoirsTest


class TestPDHandler(PynterTest):
    
    def setUp(self):
        self.comp = Composition('NaNbO3')
        self.pd = get_object_from_json(PhaseDiagram,self.get_testfile_path('PD_Na-Nb-O.json'))
        self.pdh = PDHandler(self.pd)
        
    def test_calculate_single_chempot(self):
        self.assert_all_close( self.pdh.calculate_single_chempot(self.comp, Chempots({'O':-1.92,'Na':-2.26})) , -6.26 ,atol=1e-02)
        
    def test_boundary_chempots(self):
        boundary_chempots = {
         'NaNbO3-Na2Nb3O6-Nb12O29': Chempots({'O': -8.95, 'Nb': -10.86, 'Na': -2.84}),
         'NaNbO3-Na3NbO4-O2': Chempots({'O': -4.95, 'Nb': -21.43, 'Na': -4.26}),
         'NaNbO3-Na2Nb3O6-NaNbO2': Chempots({'O': -9.08, 'Nb': -10.86, 'Na': -2.44}),
         'NaNbO3-Na3NbO4-NaNbO2': Chempots({'O': -9.08, 'Nb': -11.1, 'Na': -2.2}),
         'NaNb3O8-NaNbO3-Nb12O29': Chempots({'O': -8.14, 'Nb': -12.81, 'Na': -3.31}),
         'NaNb3O8-NaNbO3-O2': Chempots({'O': -4.95, 'Nb': -20.79, 'Na': -4.91})
         }
        ReservoirsTest().assert_Reservoirs_equal( 
            Reservoirs(self.pdh.get_all_boundaries_chempots(self.comp)) , Reservoirs(boundary_chempots), check_reference=False, rtol=1e-02 )
        
    def test_entries(self):
        comp , pd , pdh = self.comp, self.pd, self.pdh
        target_entries=[]
        for e in pd.all_entries:
            if e.composition.reduced_composition == comp:
                target_entries.append(e)
        assert pdh.get_entries_from_comp(comp) == target_entries
        assert pdh.get_stable_entry_from_comp(comp) == pd.entries[94]
        self.assert_all_close (pdh.get_formation_energy_from_stable_comp(comp), -14.28 ,atol=1e-02)
        
    def test_get_phase_boundary_chempots(self):
        phase_boundary_chempots = {
         'NaNb3O8-NaNbO3': Chempots({'Na': -2.64, 'Nb': -5.89, 'O': -1.92}),
         'NaNbO3-Na3NbO4': Chempots({'Na': -1.99, 'Nb': -6.53, 'O': -1.92})
         }
        ReservoirsTest().assert_Reservoirs_equal(  
            Reservoirs(self.pdh.get_phase_boundaries_chempots(self.comp,{'O':-1.92}) ) ,
            Reservoirs(phase_boundary_chempots),
            check_reference=False, rtol=1e-02)