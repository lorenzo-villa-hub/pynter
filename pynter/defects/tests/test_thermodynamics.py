#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:56:05 2023

@author: villa
"""
import os
import os.path as op
import unittest
from scipy import constants
from scipy.constants import physical_constants 
import numpy as np
from pymatgen.electronic_structure.dos import CompleteDos

from pynter.defects.analysis import DefectsAnalysis
from pynter.defects.thermodynamics import PressureAnalysis, ThermoData, Conductivity
from pynter.phase_diagram.chempots import PressureReservoirs
from pynter.tools.utils import get_object_from_json

from pynter.testing.core import PynterTest

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/defects/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)


e = constants.e
kb = physical_constants['Boltzmann constant in eV/K'][0]


class TestConductivity(PynterTest):
    
    def setUp(self):
        self.T = 1300
        barriers = {'Vac_Na':2.5,'Vac_O':0.8}
        mobilities = {}
        mobilities['electrons'], mobilities['holes'] = self.mobility_e(self.T), self.mobility_h(self.T)
        for name,barrier in barriers.items():
            mobilities[name] = self.mobility_ionic(self.T,barrier)
        self.conductivity = Conductivity(mobilities)
            
    def mobility_e(self,temperature):
        return 8080*(temperature**-1.5)*np.exp(-0.021/(kb*temperature))

    def mobility_h(self,temperature):
        return self.mobility_e(temperature)*0.5

    def mobility_ionic(self,temperature,barrier):
        D0 = 1e-03 # cubic crystals
        diffusivity = D0 * np.exp(-1*barrier/(kb*temperature))
        return (e*diffusivity)/(kb*temperature)   

    def test_get_conductivity(self):
        carrier_concentrations = (2.1589205669441444e+16, 1.7962030895263843e+17)
        defect_concentrations = ThermoData.from_json(get_path('thermodata_NN_1300K.json')).defect_concentrations[2]
        sigma_test = 0.004360094243796292
        sigma = self.conductivity.get_conductivity(carrier_concentrations, defect_concentrations)
        assert sigma == sigma_test
        
        
        
class TestPressureAnalysis(unittest.TestCase):
    
    def setUp(self):
        self.da = DefectsAnalysis.from_json(get_path('DA_NN_vacancies.json'))
        self.pres = PressureReservoirs.from_json(get_path('pressure_reservoirs_NN.json'))
        self.dos = get_object_from_json(CompleteDos, get_path('Si_DOS.json')) 
        self.pa = PressureAnalysis(defects_analysis=self.da, bulk_dos=self.dos)
        self.thermo_test = ThermoData.from_json(get_path('thermodata_NN_1300K.json'))
        
    def test_get_concentrations(self):
        thermo = self.pa.get_concentrations(self.pres,temperature=1300,name='NN-test')
        assert thermo.as_dict() == self.thermo_test.as_dict()
        
    def test_get_conductivities(self):
        tc = TestCoductivity()
        tc.setUp()
        mobilities = tc.conductivity.mobilities
        thermo_cond = self.pa.get_conductivities(self.pres,mobilities,1300,'NN-sigma-test')
        conductivities_test = [
         21.975733745739,
         0.3064336399157439,
         0.004360094243796292,
         0.018329168707032613,
         0.7672079995432416
         ]
        assert thermo_cond.conductivities == conductivities_test
        
    def test_get_fermi_levels(self):
        thermo = self.pa.get_fermi_levels(self.pres,1300,'NN-test')
        assert thermo.fermi_levels == self.thermo_test.fermi_levels
        
    def test_get_quenched_fermi_levels(self):
        thermoq = self.pa.get_quenched_fermi_levels(self.pres, 1300, 300)
        fermi_levels = [
         2.4381471172332763,
         1.9536618465423583,
         1.8211702342987057,
         0.0953043087005615,
         0.0018720005035399972
         ]
        assert thermoq.fermi_levels == fermi_levels
        
        thermoq = self.pa.get_quenched_fermi_levels(self.pres,1300,300,quenched_species=['Vac_Na'])
        fermi_levels = [
         2.096203606796265,
         1.5317770740509031,
         1.047121461868286,
         0.5623636447906493,
         0.07581724205017086
         ]
        assert thermoq.fermi_levels == fermi_levels
        
        thermoq = self.pa.get_quenched_fermi_levels(self.pres,1300,300,quenched_species=['Vac_O'])
        fermi_levels = [
         2.454942788314819,
         2.33268869972229,
         1.8735332088470458,
         1.3932383388519287,
         0.8850756008148193
         ]
        
        
        