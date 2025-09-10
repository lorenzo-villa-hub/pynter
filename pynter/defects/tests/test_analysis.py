#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 17:34:15 2023

@author: villa
"""
import matplotlib

#matplotlib.use('Agg') # no graphical output

from pymatgen.electronic_structure.dos import CompleteDos
import numpy as np
import os.path as op

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots
from pynter.defects.defects import DefectName, DefectComplexName
from pynter.defects.analysis import DefectsAnalysis, SingleDefConc, DefectConcentrations

from pynter.testing.core import PynterTest
from pynter.testing.defects import DefectEntryTest
from pynter.defects.tests.test_entries import TestDefectEntry


def assert_concentration_equal(actual_concentration_object,desired_charge,desired_conc,desired_name):
    sdc = actual_concentration_object
    assert sdc.charge == float(desired_charge)
    np.testing.assert_allclose(sdc.conc, desired_conc, rtol=100, atol=0.2) 
    assert sdc.name == desired_name


class TestDefectsAnalysis(PynterTest):
    
    @classmethod
    def setUpClass(cls):
        self = cls()
        da = DefectsAnalysis.from_vasp_directories(
                                                    path_defects=op.join(self.test_files_path,'Si-defects/Defects'),
                                                    path_bulk=op.join(self.test_files_path,'Si-defects/Bulk/Bulk-3x3x3-supercell'),
                                                    common_path='2-PBE-OPT',
                                                    tol=5e-03)
        da.filter_entries(inplace=True,exclude=True,types=['DefectComplex'])
        for idx in [0,2,4]:
            da[idx].label = 'mult108'
        for idx in [1,3,5]:
            da[idx].label = 'mult54'
        cls.da = da
        mu_B = -6.6794
        mu_P = -5.4133
        mu_Si = -5.4224
        cls.chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B},ndecimals=2)
        cls.dos = get_object_from_json(CompleteDos, self.get_testfile_path('Si_DOS.json')) 


    def test_import(self):
        assert self.da.occupation_function == 'MB' #check MB is default    
        entry = TestDefectEntry().get_entry()
        DefectEntryTest().assert_DefectEntry_equal(self.da.entries[3], entry)
        
    def test_stable_charges(self):
        stable_charges = {
            'Int_Si(mult108)': (0.0, 3.8291055299999766),
            'Int_Si(mult54)': (1.0, 2.909830309999954),
            'Sub_B_on_Si': (-1.0, 0.5163766299999679),
            'Sub_P_on_Si': (1.0, -0.4067508600000451),
            'Vac_Si': (0.0, 3.2633273699999723)}
        self.assert_object_almost_equal(self.da.stable_charges(self.chempots), stable_charges)
    
    def test_formation_energies(self):
        formation_energies = {
            'Int_Si(mult108)': [(-1.0, 4.700175059999967),
            (-1.0, 4.761992819999964),
            (0.0, 3.8291055299999766)],
            'Int_Si(mult54)': [(0.0, 3.8109430399999678),
            (1.0, 3.380026319999949),
            (1.0, 2.909830309999954)],
            'Sub_B_on_Si': [(-2.0, 1.127356349999948),
            (-1.0, 0.5163766299999679),
            (0.0, 0.7012506599999799)],
            'Sub_P_on_Si': [(0.0, 0.2571857899999852),
            (1.0, -0.4067508600000451),
            (2.0, -0.25286366000005067)],
            'Vac_Si': [(-1.0, 3.3212670499999675),
            (0.0, 3.2633273699999723),
            (1.0, 3.4789555099999596)]}
        
        self.assert_object_almost_equal( self.da.formation_energies(self.chempots) , formation_energies )
    
    def test_charge_transition_levels(self):
        charge_transition_levels = {
            'Int_Si(mult108)': [(0, -1, 0.8715687999999413)],
            'Int_Si(mult54)': [(1, 0, 0.90145389999994)],
            'Sub_B_on_Si': [(0, -1, -0.1843714000000135), (-1, -2, 0.6109916999999525)],
            'Sub_P_on_Si': [(2, 1, -0.15343770000001483), (1, 0, 0.6639459999999502)],
            'Vac_Si': [(1, 0, -0.21530510000001218), (0, -1, 0.05837949999997605)]}
        
        self.assert_object_almost_equal( self.da.charge_transition_levels() , charge_transition_levels )
    
    def test_carrier_concentrations(self):
        carrier_concentrations = (4.757391414506028e+19, 8150632633.6063385)
        self.assert_all_close( self.da.carrier_concentrations(self.dos) , carrier_concentrations )
    
    def test_defect_concentrations(self):
        conc =  self.da.defect_concentrations(self.chempots,fermi_level=0.4)
        assert_concentration_equal(conc[0], desired_charge=-1, desired_conc=5.33e-52, desired_name='Int_Si(mult108)')
        assert_concentration_equal(conc[1], desired_charge=-1, desired_conc=4.87e-53, desired_name='Int_Si(mult108)')
        assert_concentration_equal(conc[2], desired_charge=0, desired_conc=4.37e-44, desired_name='Int_Si(mult108)')
        assert_concentration_equal(conc[3], desired_charge=0, desired_conc=8.81e-44, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[4], desired_charge=1, desired_conc=2.91e-43, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[5], desired_charge=1, desired_conc=2.31e-35, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[6], desired_charge=-2, desired_conc=2.93e+15, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[7], desired_charge=-1, desired_conc=1.03e+19, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[8], desired_charge=0, desired_conc=1.53e+09, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[9], desired_charge=0, desired_conc=4.42e+16, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[10], desired_charge=1, desired_conc=1.20e+21, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[11], desired_charge=2, desired_conc=5.95e+11, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[12], desired_charge=-1, desired_conc=7.78e-29, desired_name='Vac_Si')
        assert_concentration_equal(conc[13], desired_charge=0, desired_conc=1.40e-34, desired_name='Vac_Si')
        assert_concentration_equal(conc[14], desired_charge=1, desired_conc=6.35e-45, desired_name='Vac_Si')

        
    def test_defect_concentrations_fixed(self):    
        concfix = 1e17
        fixed = {'P':concfix,'B':concfix}
        conc = self.da.defect_concentrations(self.chempots,fermi_level=0.4,fixed_concentrations=fixed)
        assert_concentration_equal(conc[0], desired_charge=-1, desired_conc=5.33e-52, desired_name='Int_Si')
        assert_concentration_equal(conc[1], desired_charge=-1, desired_conc=4.87e-53, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[2], desired_charge=0, desired_conc=4.37e-44, desired_name='Int_Si')
        assert_concentration_equal(conc[3], desired_charge=0, desired_conc=8.81e-44, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[4], desired_charge=1, desired_conc=2.91e-43, desired_name='Int_Si')
        assert_concentration_equal(conc[5], desired_charge=1, desired_conc=2.31e-35, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[6], desired_charge=-2, desired_conc=2.85e+13, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[7], desired_charge=-1, desired_conc=9.97e+16, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[8], desired_charge=0, desired_conc=1.49e+07, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[9], desired_charge=-1, desired_conc=1.74e+07, desired_name='Sub_B_on_Si-Sub_P_on_Si')
        assert_concentration_equal(conc[10], desired_charge=0, desired_conc=1.92e+10, desired_name='Sub_B_on_Si-Sub_P_on_Si')
        assert_concentration_equal(conc[11], desired_charge=1, desired_conc=8.33e-01, desired_name='Sub_B_on_Si-Sub_P_on_Si')
        assert_concentration_equal(conc[12], desired_charge=0, desired_conc=3.68e+12, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[13], desired_charge=1, desired_conc=1.00e+17, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[14], desired_charge=2, desired_conc=4.96e+07, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[15], desired_charge=-1, desired_conc=7.78e-29, desired_name='Vac_Si')
        assert_concentration_equal(conc[16], desired_charge=0, desired_conc=1.40e-34, desired_name='Vac_Si')
        assert_concentration_equal(conc[17], desired_charge=1, desired_conc=6.35e-45, desired_name='Vac_Si')

        
        conc_elemental = {
            'Si': 2.3092704698607755e-35,
            'B': 9.999999999999998e+16,
            'P': 1e+17,
            'Vac_Si': 7.781280016145727e-29}
        self.assert_object_almost_equal(conc.elemental , conc_elemental,rtol=100,atol=1e-02 )
        
        conc_total = {
            'Int_Si(mult108)': 4.365510595961556e-44,
            'Int_Si(mult54)': 2.309270465495265e-35,
            'Sub_B_on_Si': 1e+17,
            'Sub_P_on_Si': 1e+17,
            'Vac_Si': 7.781280016145727e-29}
        
        self.assert_object_almost_equal( conc.total , conc_total ,rtol=100,atol=1e-02)
    
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
        fermi_level = 0.46143868598937976
        self.assert_all_close( self.da.solve_fermi_level(self.chempots,self.dos) , fermi_level )

    
    def test_sort_entries(self):
        desired_indexes = [6, 0, 1, 7, 12, 2, 3, 8, 9, 13, 4, 5, 10, 14, 11]
        actual_indexes = []    
        for entry in self.da.sort_entries(features=['charge']):
            actual_indexes.append(self.da.entries.index(entry))

        self.assertEqual(actual_indexes, desired_indexes)
        
              
class TestDefectsAnalysisComplexes(PynterTest):
     
    @classmethod
    def setUpClass(cls):
        self = cls()
        da = DefectsAnalysis.from_vasp_directories(
                                                    path_defects=op.join(self.test_files_path,'Si-defects/Defects'),
                                                    path_bulk=op.join(self.test_files_path,'Si-defects/Bulk/Bulk-3x3x3-supercell'),
                                                    common_path='2-PBE-OPT',
                                                    tol=5e-03)
        for idx in [1,3,5]:
            da[idx].label = 'mult54'
        cls.da = da
        mu_B = -6.6794
        mu_P = -5.4133
        mu_Si = -5.4224
        cls.chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B},ndecimals=2)
        
    def test_binding_energy(self):
        self.assert_all_close( self.da.binding_energy('Sub_B_on_Si-Sub_P_on_Si') , 0.17095738713175915 )
        # include different labels
        
    def test_defect_concentrations_fixed(self):
        concfix = 1e17
        fixed = {'P':concfix}
        conc = self.da.defect_concentrations(self.chempots,fixed_concentrations=fixed,fermi_level=0.4)
        assert_concentration_equal(conc[0], desired_charge=-1, desired_conc=5.33e-52, desired_name='Int_Si')
        assert_concentration_equal(conc[1], desired_charge=-1, desired_conc=4.87e-53, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[2], desired_charge=0, desired_conc=4.37e-44, desired_name='Int_Si')
        assert_concentration_equal(conc[3], desired_charge=0, desired_conc=8.81e-44, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[4], desired_charge=1, desired_conc=2.91e-43, desired_name='Int_Si')
        assert_concentration_equal(conc[5], desired_charge=1, desired_conc=2.31e-35, desired_name='Int_Si(mult54)')
        assert_concentration_equal(conc[6], desired_charge=-2, desired_conc=2.85e+13, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[7], desired_charge=-1, desired_conc=9.97e+16, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[8], desired_charge=0, desired_conc=1.49e+07, desired_name='Sub_B_on_Si')
        assert_concentration_equal(conc[9], desired_charge=-1, desired_conc=1.74e+07, desired_name='Sub_B_on_Si-Sub_P_on_Si')
        assert_concentration_equal(conc[10], desired_charge=0, desired_conc=1.92e+10, desired_name='Sub_B_on_Si-Sub_P_on_Si')
        assert_concentration_equal(conc[11], desired_charge=1, desired_conc=8.33e-01, desired_name='Sub_B_on_Si-Sub_P_on_Si')
        assert_concentration_equal(conc[12], desired_charge=0, desired_conc=3.68e+12, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[13], desired_charge=1, desired_conc=1.00e+17, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[14], desired_charge=2, desired_conc=4.96e+07, desired_name='Sub_P_on_Si')
        assert_concentration_equal(conc[15], desired_charge=-1, desired_conc=7.78e-29, desired_name='Vac_Si')
        assert_concentration_equal(conc[16], desired_charge=0, desired_conc=1.40e-34, desired_name='Vac_Si')
        assert_concentration_equal(conc[17], desired_charge=1, desired_conc=6.35e-45, desired_name='Vac_Si')

        
        conc_elemental = {
            'Si': 2.3092704698607755e-35,
            'B': 9.976863341411933e+16,
            'P': 9.999803777500397e+16,
            'Vac_Si': 7.781280016145727e-29}
        self.assert_object_almost_equal( conc.elemental , conc_elemental, rtol=100,atol=1e-02 )
        
        conc_total = {
                'Int_Si': 3.9784767205219347e-44,
                'Int_Si(mult54)': 2.104537049311212e-35,
                'Sub_B_on_Si': 1.1525432095400974e+19,
                'Sub_B_on_Si-Sub_P_on_Si': 2225298056161.7256,
                'Sub_P_on_Si': 9.999777470194382e+16,
                'Vac_Si': 8.538257918013443e-29}
        self.assert_object_almost_equal( conc.total , conc_total ,rtol=100,atol=1e-02)
        
    
    def test_names(self):
        desired_names = [
            DefectName.from_string("Int_Si").as_dict(),
            DefectName.from_string("Int_Si").as_dict(),
            DefectName.from_string("Int_Si").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
            DefectName.from_string("Int_Si(mult54)").as_dict(),
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
    

    @classmethod
    def setUpClass(cls):
        self = cls()
        da = DefectsAnalysis.from_vasp_directories(
                                                    path_defects=op.join(self.test_files_path,'Si-defects/Defects'),
                                                    path_bulk=op.join(self.test_files_path,'Si-defects/Bulk/Bulk-3x3x3-supercell'),
                                                    common_path='2-PBE-OPT',
                                                    tol=5e-03)
        for idx in [3,4,5]:
            da[idx].label = 'mult54'
        cls.da = da
        mu_B = -6.6794
        mu_P = -5.4133
        mu_Si = -5.4224
        cls.chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B},ndecimals=2)
        concfix = 1e17
        fixed = {'P':concfix}
        self.conc = da.defect_concentrations(cls.chempots,fixed_concentrations=fixed,fermi_level=0.4)
        
    def test_as_dict_from_dict(self):
        conc_test = DefectConcentrations.from_dict(self.conc.as_dict())
        self.assert_object_almost_equal( self.conc.as_dict() , conc_test.as_dict() )
        
    def test_total(self):
        total = {
            'Int_Si': 3.9784767205219347e-44,
            'Int_Si(mult54)': 2.104537049311212e-35,
            'Sub_B_on_Si': 1.1525432095400974e+19,
            'Sub_B_on_Si-Sub_P_on_Si': 2225298056161.7256,
            'Sub_P_on_Si': 9.999777470194382e+16,
            'Vac_Si': 8.538257918013443e-29}
        self.assert_object_almost_equal( self.conc.total , total ,rtol=100,atol=1e-02)
        
    def test_elemental(self):
        elemental_desired = {
            'Si': 2.1045370532896886e-35,
            'B': 1.152543432069903e+19,
            'P': 9.999999999999998e+16,
            'Vac_Si': 8.538257918013443e-29}
        self.assert_object_almost_equal( self.conc.elemental, elemental_desired, rtol=100,atol=1e-02)
        
    def test_stable(self):
        stable =self.conc.stable
        assert_concentration_equal(stable[0], desired_charge=0, desired_conc=3.98e-44, desired_name='Int_Si')
        assert_concentration_equal(stable[1], desired_charge=1, desired_conc=2.10e-35, desired_name='Int_Si(mult54)')
        assert_concentration_equal(stable[2], desired_charge=-1, desired_conc=1.15e+19, desired_name='Sub_B_on_Si')
        assert_concentration_equal(stable[3], desired_charge=0, desired_conc=2.22e+12, desired_name='Sub_B_on_Si-Sub_P_on_Si')
        assert_concentration_equal(stable[4], desired_charge=1, desired_conc=1.00e+17, desired_name='Sub_P_on_Si')
        assert_concentration_equal(stable[5], desired_charge=-1, desired_conc=8.54e-29, desired_name='Vac_Si')


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
        
    
    
    
    
        
        
        
        
        
        
        
        
        
        
        
        
    