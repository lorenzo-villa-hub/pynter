#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:02:03 2023

@author: villa
"""
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots, Reservoirs, PDHandler

from pynter.testing.core import PynterTest
from pynter.testing.phase_diagram import ChempotsTest, ReservoirsTest


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
        ReservoirsTest().assert_Reservoirs_equal( self.res , 
                                        Reservoirs.from_json(self.get_testfile_path('reservoirs_NaNbO3.json')), rtol=1e-02)
    
    def test_absolute_referenced(self):
        self.res.set_to_absolute()
        ReservoirsTest().assert_Reservoirs_equal( self.res , Reservoirs({'X':self.mu_abs}) , check_reference=False, rtol=1e-03 )
        self.res.set_to_referenced()
        ReservoirsTest().assert_Reservoirs_equal( self.res , Reservoirs({'X':self.mu}), check_reference=False, rtol=1e-03)
        res_filtered = Reservoirs({'X':Chempots({'Na':-2.26,'O':-1.92})})
        ReservoirsTest().assert_Reservoirs_equal( 
            self.res.filter_reservoirs(elements=['Na','O']), res_filtered , check_reference=False, rtol=1e-03)
    
    
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
    
    
    
    
    
    
    
    
    
    