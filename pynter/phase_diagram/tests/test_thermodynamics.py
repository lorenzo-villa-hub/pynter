#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 11:11:39 2023

@author: villa
"""
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core.composition import Composition

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots, PressureReservoirs
from pynter.phase_diagram.thermodynamics import OxygenPressure

from pynter.testing.phase_diagram import ReservoirsTest
from pynter.phase_diagram.tests.test_chempots import TestChempots


class TestOxygenPressure(TestChempots):
    
    def setUp(self):
        super().setUp()
        self.oxpr = OxygenPressure(temperature=1300,partial_pressure=0.2)
        
        
    def test_chempot_ideal_gas(self):
        oxygen_standard_chempot = self.oxpr.oxygen_standard_chempot()
        self.assert_all_close( oxygen_standard_chempot , -1.4265454545454541 )
        self.assert_all_close( self.oxpr.chempot_ideal_gas(oxygen_standard_chempot,1300) , -1.516694354091398 )
        
    def test_chempot_from_pO2(self):
        self.assert_all_close( self.oxpr.get_oxygen_chempot_from_pO2() , -1.516694354091398 )
        
    def test_oxygen_pressure_reservoirs(self):
        ox_res_dict = {
          1e-20: Chempots({'O': -10.88}),
          3.162277660168379e-13: Chempots({'O': -9.91}),
          1e-05: Chempots({'O': -8.94}),
          316.22776601683796: Chempots({'O': -7.97}),
          10000000000.0: Chempots({'O': -7.01})
          }
        ox_res = self.oxpr.get_oxygen_pressure_reservoirs(self.mu_abs['O'],npoints=5)
        print(ox_res)
        ReservoirsTest().assert_Reservoirs_equal( ox_res , PressureReservoirs(ox_res_dict) ,check_reference=False)
        
    def test_pressure_reservoirs_from_pd(self):
        pres_dict = {
         1e-20: {'Na': -2.53, 'Nb': -11.14, 'O': -8.96},
         3.162277660168379e-13: {'Na': -3.07, 'Nb': -13.51, 'O': -7.99},
         1e-05: {'Na': -3.55, 'Nb': -15.94, 'O': -7.02},
         316.22776601683796: {'Na': -4.03, 'Nb': -18.36, 'O': -6.05},
         10000000000.0: {'Na': -4.51, 'Nb': -20.76, 'O': -5.09}
         }
        pd = get_object_from_json(PhaseDiagram,self.get_testfile_path('PD_Na-Nb-O.json'))
        comp = Composition('NaNbO3')
        pres = self.oxpr.get_pressure_reservoirs_from_pd(pd,comp,temperature=1300,npoints=5)
        ReservoirsTest().assert_Reservoirs_equal( pres , PressureReservoirs(pres_dict) ,check_reference=False)
        ReservoirsTest().assert_Reservoirs_equal( pres , PressureReservoirs.from_dict(pres.as_dict()) ,check_reference=False)
    