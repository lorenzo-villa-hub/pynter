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
from pynter.defects.chempots.core import Chempots
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
                                                    path_defects=op.join(self.test_files_path,'SiO2-defects/Defects'),
                                                    path_bulk=op.join(self.test_files_path,'SiO2-defects/Bulk-2x2x2-supercell'),
                                                    get_charge_correction=False,
                                                    common_path='2-PBE-OPT',
                                                    initial_structure=True)
        cls.da_comp = da.copy()
        da.filter_entries(inplace=True,exclude=True,types=['DefectComplex'])
        cls.da = da
        mu_O = -6
        mu_P = -5.4133
        mu_Si = -5.4224
        cls.chempots = Chempots({'Si':mu_Si,'P':mu_P,'O':mu_O},ndecimals=2)
        cls.dos = get_object_from_json(CompleteDos, self.get_testfile_path('SiO2-defects/Bulk-2x2x2-supercell/complete_dos.json')) 


        
    def test_stable_charges(self):
        stable_charge_Vac_O = (2, 0.08920299999999681)
        stable_charge_Int_O = (0, 5.666554999999903)
        self.assert_all_close(self.da.stable_charges(self.chempots)['Vac_O'], stable_charge_Vac_O)
        self.assert_all_close(self.da.stable_charges(self.chempots)['Int_O'], stable_charge_Int_O)
    
    def test_formation_energies(self):
        formation_energies = self.da.formation_energies(self.chempots)
        formation_energy_Vac_Si_qminus4 = 21.89906667999999
        desired = formation_energy_Vac_Si_qminus4
        actual = formation_energies['Vac_Si'][1][1]
        self.assert_all_close(actual , desired, rtol=1e-03 )
    
    def test_charge_transition_levels(self):
        actual = self.da.charge_transition_levels()['Sub_P_on_Si'][0]
        desired = (1, 0, 6.224602499999989)
        self.assert_all_close( actual , desired )
    
    def test_carrier_concentrations(self):
        actual = self.da.carrier_concentrations(bulk_dos=self.dos,fermi_level=2)
        desired = (2.5207523856086593e-16, 1.9354062962348537e-41)
        self.assert_all_close(actual , desired, rtol=1e-03 )
    
    def test_defect_concentrations(self):
        conc = self.da.defect_concentrations(self.chempots)
        actual = conc.select_concentrations(name='Vac_O',charge=2)[0]['conc']
        desired = 3.383174303005728e+19
        self.assert_all_close(actual, desired, rtol=1e-03)
        actual = conc.select_concentrations(name='Sub_P_on_Si',charge=1)[0]['conc']
        desired = 1.1001057398319168e+21
        self.assert_all_close(actual, desired, rtol=1e-03)

        
    def test_defect_concentrations_fixed(self):    
        concfix = 1e17
        fixed = {'P':concfix}
        conc = self.da.defect_concentrations(chemical_potentials=self.chempots,fixed_concentrations=fixed)

        actual = conc.elemental['P']
        desired = concfix
        self.assert_all_close(actual, desired)

        actual = conc.select_concentrations(charge=2, name='Vac_O')[0]['conc']
        desired = 3.383174303005728e+19
        self.assert_all_close(actual, desired, rtol=1e-03)

        fixed = {'Vac_O':concfix}
        conc = self.da.defect_concentrations(chemical_potentials=self.chempots,fixed_concentrations=fixed,fermi_level=1.87)
        conc_Vac_O_q2 = conc.select_concentrations(name='Vac_O',charge=2)[0]['conc']
        conc_Vac_O_q0 = conc.select_concentrations(name='Vac_O',charge=0)[0]['conc']
        actual = conc_Vac_O_q2, conc_Vac_O_q0
        desired = (6.412340108578859e+16,3.58765989142114e+16)
        self.assert_all_close(actual, desired, rtol=1e-03)
        self.assert_all_close(conc.elemental['Vac_O'], concfix, rtol=1e-05)
        self.assert_all_close(conc.total['Sub_P_on_Si'], 0.13903596778575308, rtol=5e-03)

    def test_names(self):
        actual = self.da.names
        desired = ['Int_O', 'Sub_P_on_Si', 'Vac_O', 'Vac_Si', 'Vac_Si-Vac_O']
        self.assert_object_almost_equal(actual, desired)
        
    def test_select_entries(self):
        manual_entries = self.da[3], self.da[4]
        selected_entries = self.da.select_entries(elements=['P'])
        for i in range(0,len(selected_entries)):
            DefectEntryTest().assert_DefectEntry_equal(selected_entries[i], manual_entries[i])
            
    def test_solve_fermi_level(self):
        actual = self.da.solve_fermi_level(chemical_potentials=self.chempots,bulk_dos=self.dos)
        desired = 3.0042287011146547
        self.assert_all_close(actual,desired)

    def test_binding_energy(self):
        actual = self.da_comp.binding_energy(name='Vac_Si-Vac_O')
        desired = 1.3470440099999914
        self.assert_all_close(actual, desired)
           
        
    def test_from_dataframe(self):
        da = DefectsAnalysis.from_file(
                                self.get_testfile_path('SiO2_defects_analysis.csv'),
                                vbm=self.da.vbm,
                                band_gap=self.da.band_gap)
        self.da = da
        self.test_defect_concentrations_fixed()
        self.test_formation_energies()
        self.test_charge_transition_levels()


        
    