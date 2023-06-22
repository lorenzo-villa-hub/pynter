#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 17:34:15 2023

@author: villa
"""
import matplotlib

matplotlib.use('Agg') # no graphical output

from pymatgen.electronic_structure.dos import CompleteDos

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots
from pynter.defects.defects import DefectName, DefectComplexName
from pynter.defects.analysis import DefectsAnalysis, SingleDefConc, DefectConcentrations

from pynter.testing.core import PynterTest
from pynter.testing.defects import DefectEntryTest
from pynter.defects.tests.test_entries import TestDefectEntry



class TestDefectsAnalysis(PynterTest):
    
    def setUp(self):
        self.da = DefectsAnalysis.from_json(self.get_testfile_path('DA_Si_single.json'))
        mu_B = -6.6794
        mu_P = -5.4133
        mu_Si = -5.4224
        self.chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B},round_values=2)
        self.dos = get_object_from_json(CompleteDos, self.get_testfile_path('Si_DOS.json')) 

    def test_import(self):
        assert self.da.occupation_function == 'MB' #check MB is default    
        entry = TestDefectEntry().get_entry()
        DefectEntryTest().assert_DefectEntry_equal(self.da.entries[3], entry)
        
    def test_stable_charges(self):
        stable_charges = {
             "Int_Si(mult108)": (1.0, 3.6846594634256036),
             "Int_Si(mult54)": (1.0, 3.2615360171135395),
             "Sub_B_on_Si": (0.0, 0.7012506599999799),
             "Sub_P_on_Si": (1.0, -0.25706887286827396),
             "Vac_Si": (0.0, 3.2633273699999723)
             }
        self.assert_object_almost_equal(self.da.stable_charges(self.chempots), stable_charges)
    
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
        
        self.assert_object_almost_equal( self.da.formation_energies(self.chempots) , formation_energies )
    
    def test_charge_transition_levels(self):
        charge_transition_levels = {
            'Sub_P_on_Si': [(2, 1, -0.44757000000000224), (1, 0, 0.5145204999999566)],
            'Vac_Si': [(1, 0, -0.08527870000001775), (0, -1, 0.43482689999996005)],
            'Sub_B_on_Si': [(0, -1, 0.09193469999997461)],
            'Int_Si(mult108)': [(1, 0, 0.1448889999999724), (0, -1, 0.8207116999999435)],
            'Int_Si(mult54)': [(1, 0, 0.549648599999955), (0, -1, 0.8631799999999417)]
            }
        
        self.assert_object_almost_equal( self.da.charge_transition_levels() , charge_transition_levels )
    
    def test_carrier_concentrations(self):
        carrier_concentrations = (4.757391414506028e+19, 8150632633.6063385)
        self.assert_all_close( self.da.carrier_concentrations(self.dos) , carrier_concentrations )
    
    def test_defect_concentrations(self):
        concentrations = DefectConcentrations([
    SingleDefConc(charge=-1.0, conc=4.04e-49, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
    SingleDefConc(charge=0.0, conc=4.71e-42, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
    SingleDefConc(charge=1.0, conc=2.40e-46, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
    SingleDefConc(charge=-1.0, conc=8.00e-50, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
    SingleDefConc(charge=0.0, conc=4.76e-42, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
    SingleDefConc(charge=1.0, conc=1.54e-39, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
    SingleDefConc(charge=-2.0, conc=1.99e+04, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
    SingleDefConc(charge=-1.0, conc=1.24e+16, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
    SingleDefConc(charge=0.0, conc=8.28e+10, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
    SingleDefConc(charge=0.0, conc=2.39e+18, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
    SingleDefConc(charge=1.0, conc=1.98e+20, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
    SingleDefConc(charge=2.0, conc=1.14e+06, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
    SingleDefConc(charge=-1.0, conc=1.97e-33, name=DefectName.from_string("Vac_Si"), stable=True),
    SingleDefConc(charge=0.0, conc=7.54e-33, name=DefectName.from_string("Vac_Si"), stable=True),
    SingleDefConc(charge=1.0, conc=5.20e-41, name=DefectName.from_string("Vac_Si"), stable=True)
        ])
        actual_concentrations =  self.da.defect_concentrations(self.chempots,fermi_level=0.4)

        self.assert_object_almost_equal(actual_concentrations.as_dict() ,
                                      concentrations.as_dict(),rtol=100,atol=1e-02) #ignore rtol
    
    
    def test_defect_concentrations_fixed(self):    
        concfix = 1e17
        fixed = {'P':concfix,'B':concfix}
        conc_fixed = self.da.defect_concentrations(self.chempots,fermi_level=0.4,fixed_concentrations=fixed)
        
        conc_fixed = DefectConcentrations([
            SingleDefConc(charge=-1.0, conc=4.04e-49, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
            SingleDefConc(charge=0.0, conc=4.71e-42, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
            SingleDefConc(charge=1.0, conc=2.40e-46, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
            SingleDefConc(charge=-1.0, conc=8.00e-50, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
            SingleDefConc(charge=0.0, conc=4.76e-42, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
            SingleDefConc(charge=1.0, conc=1.54e-39, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
            SingleDefConc(charge=-2.0, conc=1.60e+05, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
            SingleDefConc(charge=-1.0, conc=1.00e+17, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
            SingleDefConc(charge=0.0, conc=6.67e+11, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
            SingleDefConc(charge=0.0, conc=1.19e+15, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=1.0, conc=9.88e+16, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=2.0, conc=5.69e+02, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=-1.0, conc=1.97e-33, name=DefectName.from_string("Vac_Si"), stable=True),
            SingleDefConc(charge=0.0, conc=7.54e-33, name=DefectName.from_string("Vac_Si"), stable=True),
            SingleDefConc(charge=1.0, conc=5.20e-41, name=DefectName.from_string("Vac_Si"), stable=True)
            ])

        self.assert_object_almost_equal( conc_fixed.as_dict() , conc_fixed.as_dict() ,rtol=100,atol=1e-02)
        
        conc_elemental = {
         'Si': 9.504324340294853e-33,
         'B': 1.0000000000000002e+17,
         'P': 1.0000000000000002e+17
         }
        self.assert_object_almost_equal(conc_fixed.elemental , conc_elemental,rtol=100,atol=1e-02 )
        
        conc_total = {
         'Int_Si(mult108)': 4.7149918453710235e-42,
         'Int_Si(mult54)': 1.544589116486735e-39,
         'Sub_B_on_Si': 1.0000000000000002e+17,
         'Sub_P_on_Si': 1.0000000000000002e+17,
         'Vac_Si': 9.504322790990746e-33
         }
        
        self.assert_object_almost_equal( conc_fixed.total , conc_total ,rtol=100,atol=1e-02)
    
    def test_names(self):
        desired_names = [
            DefectName.from_string("Int_Si(mult108)").as_dict(),
            DefectName.from_string("Int_Si(mult108)").as_dict(),
            DefectName.from_string("Int_Si(mult108)").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
            DefectName.from_string("Sub_B_on_Si").as_dict(),
            DefectName.from_string("Sub_B_on_Si").as_dict(),
            DefectName.from_string("Sub_B_on_Si").as_dict(),
            DefectName.from_string("Sub_P_on_Si").as_dict(),
            DefectName.from_string("Sub_P_on_Si").as_dict(),
            DefectName.from_string("Sub_P_on_Si").as_dict(),
            DefectName.from_string("Vac_Si").as_dict(),
            DefectName.from_string("Vac_Si").as_dict(),
            DefectName.from_string("Vac_Si").as_dict()
            ]
        
        actual_names = [e.name.as_dict() for e in self.da]     
        self.assert_object_almost_equal(actual_names, desired_names)
        
    def test_select_entries(self):
        manual_entries = [self.da.entries[i] for i in range(6,12)]
        selected_entries = self.da.select_entries(elements=['B','P'])
        for i in range(0,len(selected_entries)):
            DefectEntryTest().assert_DefectEntry_equal(selected_entries[i], manual_entries[i])
            
    def test_solve_fermi_level(self):
        fermi_level = 0.5058802570343016
        self.assert_all_close( self.da.solve_fermi_level(self.chempots,self.dos) , fermi_level )

    
    def test_sort_entries(self):
        indexes = [6, 0, 3, 7, 12, 1, 4, 8, 9, 13, 2, 5, 10, 14, 11]
        sorted_entries = [self.da.entries[i] for i in indexes]
        assert self.da.sort_entries(features=['charge']) == sorted_entries
        
        
        
class TestDefectsAnalysisComplexes(PynterTest):
    
    def setUp(self):
        self.da = DefectsAnalysis.from_json(self.get_testfile_path('DA_Si.json'))
        mu_B = -6.6794
        mu_P = -5.4133
        mu_Si = -5.4224
        self.chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B})
        self.dos = get_object_from_json(CompleteDos, self.get_testfile_path('Si_DOS.json'))         
        
    def test_binding_energy(self):
        self.assert_all_close( self.da.binding_energy('Int_Si-Vac_Si') , -1.5579998134255888 )
        self.assert_all_close( self.da.binding_energy('Sub_B_on_Si-Sub_P_on_Si') , -0.17095738713175915 )
        
    def test_defect_concentrations_fixed(self):
        concfix = 1e17
        fixed = {'P':concfix}
        conc = self.da.defect_concentrations(self.chempots,fixed_concentrations=fixed,fermi_level=0.4)
        
        desired_conc = DefectConcentrations([
            SingleDefConc(charge=-1.0, conc=3.68e-49, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
            SingleDefConc(charge=0.0, conc=4.30e-42, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
            SingleDefConc(charge=1.0, conc=2.19e-46, name=DefectName.from_string("Int_Si(mult108)"), stable=True),
            SingleDefConc(charge=-1.0, conc=7.29e-50, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
            SingleDefConc(charge=0.0, conc=4.34e-42, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
            SingleDefConc(charge=1.0, conc=1.40e-39, name=DefectName.from_string("Int_Si(mult54)"), stable=True),
            SingleDefConc(charge=-1.0, conc=3.96e-74, name=DefectComplexName.from_string("Int_Si-Vac_Si"), stable=True),
            SingleDefConc(charge=0.0, conc=6.56e-73, name=DefectComplexName.from_string("Int_Si-Vac_Si"), stable=True),
            SingleDefConc(charge=1.0, conc=1.08e-74, name=DefectComplexName.from_string("Int_Si-Vac_Si"), stable=True),
            SingleDefConc(charge=-2.0, conc=2.23e+04, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
            SingleDefConc(charge=-1.0, conc=1.40e+16, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
            SingleDefConc(charge=0.0, conc=9.30e+10, name=DefectName.from_string("Sub_B_on_Si"), stable=True),
            SingleDefConc(charge=-1.0, conc=2.53e+12, name=DefectComplexName.from_string("Sub_B_on_Si-Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=0.0, conc=2.79e+15, name=DefectComplexName.from_string("Sub_B_on_Si-Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=1.0, conc=1.21e+05, name=DefectComplexName.from_string("Sub_B_on_Si-Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=0.0, conc=1.16e+15, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=1.0, conc=9.60e+16, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=2.0, conc=5.53e+02, name=DefectName.from_string("Sub_P_on_Si"), stable=True),
            SingleDefConc(charge=-1.0, conc=2.16e-33, name=DefectName.from_string("Vac_Si"), stable=True),
            SingleDefConc(charge=0.0, conc=8.27e-33, name=DefectName.from_string("Vac_Si"), stable=True),
            SingleDefConc(charge=1.0, conc=5.71e-41, name=DefectName.from_string("Vac_Si"), stable=True)
        ])
        print(conc.as_dict(),'\n \n',desired_conc.as_dict())
        self.assert_object_almost_equal( conc.as_dict() , desired_conc.as_dict(),rtol=100,atol=1e-02 ,verbose=True)
        
        conc_elemental = {'Si': 9.504324340294853e-33, 'B': 1.4920566601530766e+16, 'P': 1e+17}
        self.assert_object_almost_equal( conc.elemental , conc_elemental, rtol=100,atol=1e-02 )
        
        conc_total = {
         'Int_Si(mult108)': 4.7149918453710235e-42,
         'Int_Si(mult54)': 1.544589116486735e-39,
         'Int_Si-Vac_Si': 7.062810210333003e-73,
         'Sub_B_on_Si': 1.2423427771347858e+16,
         'Sub_B_on_Si-Sub_P_on_Si': 2497138830182908.5,
         'Sub_P_on_Si': 9.750286116981709e+16,
         'Vac_Si': 9.504322790990746e-33
         }
        self.assert_object_almost_equal( conc.total , conc_total ,rtol=100,atol=1e-02)
        
    
    def test_names(self):
        desired_names = [
            DefectName.from_string("Int_Si(mult108)").as_dict(),
            DefectName.from_string("Int_Si(mult108)").as_dict(),
            DefectName.from_string("Int_Si(mult108)").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
            DefectComplexName.from_string("Int_Si-Vac_Si").as_dict(),
            DefectComplexName.from_string("Int_Si-Vac_Si").as_dict(),
            DefectComplexName.from_string("Int_Si-Vac_Si").as_dict(),
            DefectName.from_string("Sub_B_on_Si").as_dict(),
            DefectName.from_string("Sub_B_on_Si").as_dict(),
            DefectName.from_string("Sub_B_on_Si").as_dict(),
            DefectComplexName.from_string("Sub_B_on_Si-Sub_P_on_Si").as_dict(),
            DefectComplexName.from_string("Sub_B_on_Si-Sub_P_on_Si").as_dict(),
            DefectComplexName.from_string("Sub_B_on_Si-Sub_P_on_Si").as_dict(),
            DefectName.from_string("Sub_P_on_Si").as_dict(),
            DefectName.from_string("Sub_P_on_Si").as_dict(),
            DefectName.from_string("Sub_P_on_Si").as_dict(),
            DefectName.from_string("Vac_Si").as_dict(),
            DefectName.from_string("Vac_Si").as_dict(),
            DefectName.from_string("Vac_Si").as_dict()
                ]
        
        actual_names = [e.name.as_dict() for e in self.da]
        self.assert_object_almost_equal(actual_names, desired_names)
    
    def test_plot(self):
        self.da.plot(chemical_potentials=self.chempots)
        
    def test_plot_ctl(self):
        self.da.plot_ctl()
        
    def test_plot_binding_energies(self):
        self.da.plot_binding_energies()
        

class TestSingleDefConc(PynterTest):

    def setUp(self):
        self.sc = SingleDefConc(name='Vac_Si',charge=1,conc=1e10,stable=True)     

    def test_string(self):
        assert self.sc.__str__() == 'charge=1.0, conc=1.00e+10, name=Vac_Si, stable=True'
        
    def test_as_dict_from_dict(self):
        self.assertDictEqual( self.sc.as_dict() , SingleDefConc.from_dict(self.sc.as_dict()).as_dict() )
    

class TestDefectConcentrations(PynterTest):
    
    def setUp(self):
        da = DefectsAnalysis.from_json(self.get_testfile_path('DA_Si.json'))
        mu_B = -6.6794
        mu_P = -5.4133
        mu_Si = -5.4224
        chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B})
        concfix = 1e17
        fixed = {'P':concfix}
        self.conc = da.defect_concentrations(chempots,fixed_concentrations=fixed,fermi_level=0.4)
        
    def test_as_dict_from_dict(self):
        conc_test = DefectConcentrations.from_dict(self.conc.as_dict())
        self.assert_object_almost_equal( self.conc.as_dict() , conc_test.as_dict() )
        
    def test_total(self):
        total = {
         'Int_Si(mult108)': 4.7149918453710235e-42,
         'Int_Si(mult54)': 1.544589116486735e-39,
         'Int_Si-Vac_Si': 7.062810210333003e-73,
         'Sub_B_on_Si': 1.2423427771347858e+16,
         'Sub_B_on_Si-Sub_P_on_Si': 2497138830182908.5,
         'Sub_P_on_Si': 9.750286116981709e+16,
         'Vac_Si': 9.504322790990746e-33
         }
        self.assert_object_almost_equal( self.conc.total , total ,rtol=100,atol=1e-02)
        
    def test_elemental(self):
        self.assert_object_almost_equal( self.conc.elemental ,
                                        {'Si': 9.504324340294853e-33, 'B': 1.4920566601530766e+16, 'P': 1e+17} ,rtol=100,atol=1e-02)
        
    def test_stable(self):
        stable_string = (
          '[charge=0.0, conc=4.71e-42, name=Int_Si(mult108), stable=True, '
          'charge=1.0, conc=1.54e-39, name=Int_Si(mult54), stable=True, '
          'charge=0.0, conc=6.56e-73, name=Int_Si-Vac_Si, stable=True, '
          'charge=-1.0, conc=1.24e+16, name=Sub_B_on_Si, stable=True, '
          'charge=0.0, conc=2.49e+15, name=Sub_B_on_Si-Sub_P_on_Si, stable=True, '
          'charge=1.0, conc=9.63e+16, name=Sub_P_on_Si, stable=True, '
          'charge=0.0, conc=7.54e-33, name=Vac_Si, stable=True]'
          )
        self.assert_str_content_equal( self.conc.stable.__str__() , stable_string )
        
    def test_select_concentrations(self):
        sel_conc = self.conc.select_concentrations(charge=2)
        assert sel_conc[0] == self.conc[-4]
        
        sel_conc = self.conc.select_concentrations(charges=[0,1])
        test_sel_conc = [self.conc[i] for i in [1, 2, 4, 5, 7, 8, 11, 13, 14, 15, 16, 19, 20]]
        assert sel_conc == test_sel_conc
        
        sel_conc_exclude = self.conc.select_concentrations(
            sel_conc,exclude=True,mode='or',charges=[1],names=['Vac_Si','Sub_P_on_Si'])
        test_sel_conc_exclude = [self.conc[i] for i in [1, 4, 7, 11, 13]]
        assert sel_conc_exclude == test_sel_conc_exclude
        
    
    
    
    
        
        
        
        
        
        
        
        
        
        
        
        
    