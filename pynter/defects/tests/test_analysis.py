#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 17:34:15 2023

@author: villa
"""

import os
import os.path as op
import unittest
from pymatgen.electronic_structure.dos import CompleteDos

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots
from pynter.defects.analysis import DefectsAnalysis

from pynter.defects.tests.test_entries import get_defect_entry
from pynter.defects.tests.compare import CompareEntries

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/defects/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

mu_B = -6.6794
mu_P = -5.4133
mu_Si = -5.4224
chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B})
vbm = 5.8268
band_gap = 0.5729



class TestDefectsAnalysis(unittest.TestCase):
    
    def setUp(self):
        self.da = DefectsAnalysis.from_json(get_path('DA_Si_single.json'))
        mu_B = -6.6794
        mu_P = -5.4133
        mu_Si = -5.4224
        self.chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B})
        self.dos = get_object_from_json(CompleteDos, get_path('Si_DOS.json')) 

    def test_import(self):
        assert self.da.occupation_function == 'MB' #check MB is default    
        entry = get_defect_entry()
        CompareEntries().compare(self.da.entries[3], entry)
        
    def test_stable_charges(self):
        stable_charges = {
             "Int_Si(mult108)": (1.0, 3.6846594634256036),
             "Int_Si(mult54)": (1.0, 3.2615360171135395),
             "Sub_B_on_Si": (0.0, 0.7012506599999799),
             "Sub_P_on_Si": (1.0, -0.25706887286827396),
             "Vac_Si": (0.0, 3.2633273699999723)
             }
        assert self.da.stable_charges(chempots) == stable_charges 
    
    def test_formation_energies(self):
        formation_energies = {
            'Int_Si(mult108)': [(-1.0, 4.649770703425618),
             (0.0, 3.8291055299999766),
             (1.0, 3.6846594634256036)],
            'Int_Si(mult54)': [(-1.0, 4.673740193780213),
             (0.0, 3.8109430399999678),
             (1.0, 3.2615360171135395)],
            'Sub_B_on_Si': [(-2.0, 1.8952842985270144),
             (-1.0, 0.793132301342261),
             (0.0, 0.7012506599999799)],
            'Sub_P_on_Si': [(0.0, 0.2571857899999852),
             (1.0, -0.25706887286827396),
             (2.0, 0.19056955168492085)],
            'Vac_Si': [(-1.0, 3.698017458184358),
             (0.0, 3.2633273699999723),
             (1.0, 3.3491269708159415)]
             }
        
        assert self.da.formation_energies(chempots) == formation_energies
    
    def test_charge_transition_levels(self):
        charge_transition_levels = {
            'Sub_P_on_Si': [(2, 1, -0.44757000000000224), (1, 0, 0.5145204999999566)],
            'Vac_Si': [(1, 0, -0.08527870000001775), (0, -1, 0.43482689999996005)],
            'Sub_B_on_Si': [(0, -1, 0.09193469999997461)],
            'Int_Si(mult108)': [(1, 0, 0.1448889999999724), (0, -1, 0.8207116999999435)],
            'Int_Si(mult54)': [(1, 0, 0.549648599999955), (0, -1, 0.8631799999999417)]
            }
        
        assert self.da.charge_transition_levels() == charge_transition_levels
    
    def test_carrier_concentrations(self):
        carrier_concentrations = (4.757391414506028e+19, 8150632633.6063385)
        assert self.da.carrier_concentrations(self.dos) == carrier_concentrations
    
    def test_defect_concentrations(self):
        concentrations_string = (
              '[charge=-1.0, conc=4.04e-49, name=Int_Si(mult108), stable=True\n'
              'charge=0.0, conc=4.71e-42, name=Int_Si(mult108), stable=True\n'
              'charge=1.0, conc=2.40e-46, name=Int_Si(mult108), stable=True\n'
              'charge=-1.0, conc=8.00e-50, name=Int_Si(mult54), stable=True\n'
              'charge=0.0, conc=4.76e-42, name=Int_Si(mult54), stable=True\n'
              'charge=1.0, conc=1.54e-39, name=Int_Si(mult54), stable=True\n'
              'charge=-2.0, conc=1.99e+04, name=Sub_B_on_Si, stable=True\n'
              'charge=-1.0, conc=1.24e+16, name=Sub_B_on_Si, stable=True\n'
              'charge=0.0, conc=8.28e+10, name=Sub_B_on_Si, stable=True\n'
              'charge=0.0, conc=2.39e+18, name=Sub_P_on_Si, stable=True\n'
              'charge=1.0, conc=1.98e+20, name=Sub_P_on_Si, stable=True\n'
              'charge=2.0, conc=1.14e+06, name=Sub_P_on_Si, stable=True\n'
              'charge=-1.0, conc=1.97e-33, name=Vac_Si, stable=True\n'
              'charge=0.0, conc=7.54e-33, name=Vac_Si, stable=True\n'
              'charge=1.0, conc=5.20e-41, name=Vac_Si, stable=True]'
              )
        
        assert self.da.defect_concentrations(chempots,fermi_level=0.4).__str__() == concentrations_string
    
    
    def test_defect_concentrations_fixed(self):    
        site_concentration = self.da.entries[0].defect.site_concentration_in_cm3
        concfix = site_concentration*1e-05
        fixed = {'P':concfix,'B':concfix}
        
        conc_fixed_string = (
              '[charge=-1.0, conc=4.04e-49, name=Int_Si(mult108), stable=True\n'
              'charge=0.0, conc=4.71e-42, name=Int_Si(mult108), stable=True\n'
              'charge=1.0, conc=2.40e-46, name=Int_Si(mult108), stable=True\n'
              'charge=-1.0, conc=8.00e-50, name=Int_Si(mult54), stable=True\n'
              'charge=0.0, conc=4.76e-42, name=Int_Si(mult54), stable=True\n'
              'charge=1.0, conc=1.54e-39, name=Int_Si(mult54), stable=True\n'
              'charge=-2.0, conc=8.00e+05, name=Sub_B_on_Si, stable=True\n'
              'charge=-1.0, conc=5.00e+17, name=Sub_B_on_Si, stable=True\n'
              'charge=0.0, conc=3.33e+12, name=Sub_B_on_Si, stable=True\n'
              'charge=0.0, conc=5.94e+15, name=Sub_P_on_Si, stable=True\n'
              'charge=1.0, conc=4.94e+17, name=Sub_P_on_Si, stable=True\n'
              'charge=2.0, conc=2.84e+03, name=Sub_P_on_Si, stable=True\n'
              'charge=-1.0, conc=1.97e-33, name=Vac_Si, stable=True\n'
              'charge=0.0, conc=7.54e-33, name=Vac_Si, stable=True\n'
              'charge=1.0, conc=5.20e-41, name=Vac_Si, stable=True]'
              )
        
        self.da.defect_concentrations(chempots,fermi_level=0.4,fixed_concentrations=fixed) == conc_fixed_string
    
    
    def test_select_entries(self):
        manual_entries = [self.da.entries[i] for i in range(6,12)]
        selected_entries = self.da.select_entries(elements=['B','P'])
        for i in range(0,len(selected_entries)):
            CompareEntries().compare(selected_entries[i], manual_entries[i])
            
    def test_solve_fermi_level(self):
        fermi_level = 0.5058802570343016
        assert self.da.solve_fermi_level(self.chempots,self.dos) == fermi_level

    
    def test_sort_entries(self):
        indexes = [6, 0, 3, 7, 12, 1, 4, 8, 9, 13, 2, 5, 10, 14, 11]
        sorted_entries = [self.da.entries[i] for i in indexes]
        assert self.da.sort_entries(features=['charge']) == sorted_entries
        
        
        
        
        
        
        
        
        
        
        
    