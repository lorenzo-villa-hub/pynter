#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 10:15:10 2025

@author: lorenzo
"""
import numpy as np

from pynter.defects.chempots.reservoirs import PressureReservoirs
from pynter.defects.chempots.oxygen import (
                                            get_barycenter_chemical_potentials_absolute,
                                            get_barycenter_chemical_potentials_relative,
                                            get_oxygen_chempot_from_pO2,
                                            get_oxygen_chempot_standard_finite_temperature,
                                            get_oxygen_pressure_reservoirs,
                                            get_pressure_reservoirs_from_phase_diagram
                                            )
from pynter.defects.chempots.core import Chempots
from pynter.testing.chempots import ReservoirsTest

from .test_chempots import TestChempots

def test_chempot_from_pO2():
    actual = get_oxygen_chempot_from_pO2(temperature=1300,partial_pressure=0.2)
    desired = -1.516694354091398
    np.testing.assert_allclose(actual, desired, atol=1e-04)
    
    muO_reference = -5
    actual = get_oxygen_chempot_from_pO2(temperature=1300,partial_pressure=0.2,muO_reference=muO_reference)
    desired += muO_reference
    np.testing.assert_allclose(actual, desired, atol=1e-04)
    
    
    
class TestOxygenReservoirs(TestChempots):
    
    def test_oxygen_pressure_reservoirs(self):
        ox_res_dict = {
          1e-20: Chempots({'O': -10.88}),
          3.162e-13: Chempots({'O': -9.91}),
          1e-05: Chempots({'O': -8.94}),
          3.162e02: Chempots({'O': -7.97}),
          1e10: Chempots({'O': -7.01})
          }
        ox_res = get_oxygen_pressure_reservoirs(oxygen_ref=self.mu_abs['O'],temperature=1300,npoints=5)
        ReservoirsTest().assert_Reservoirs_equal( ox_res , PressureReservoirs(ox_res_dict) ,check_reference=False, rtol=1e-02)