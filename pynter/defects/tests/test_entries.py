#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 16:47:57 2023

@author: villa
"""
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots
from pynter.defects.defects import Interstitial
from pynter.defects.entries import DefectEntry

from pynter.testing.core import PynterTest


class TestDefectEntry(PynterTest):

    def setUp(self):
        mu_B = -6.6794
        mu_P = -5.4133
        mu_Si = -5.4224
        self.chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B})
        self.vbm = 5.8268
        self.band_gap = 0.5729
        self.stress_bulk = [[19.5106052, -0.0, -0.0], [0.0, 19.51060527, 0.0], [0.0, 0.0, 19.51060528]]
        self.bulk_modulus = 144
        
    @property
    def entry(self):
        return self.get_entry()
    
    def get_entry(self):
        bulk_structure = get_object_from_json(Structure,self.get_testfile_path('Si-bulk_structure_3x3x3_supercell.json'))
        defect_site = PeriodicSite('Si',[0.5417, 0.5417, 0.5417], bulk_structure.lattice)
        defect = Interstitial(defect_site, bulk_structure,charge=-1.0,multiplicity=54)
        energy_diff = 5.168792819999965
        corrections = {'kumagai': -0.08825262621975172}
        data = {'stress': 
                [[39.5772376, -0.0, 0.0],
                 [0.0, 39.57723773, 0.0],
                 [0.0, -0.0, 39.57723775]]}
        label = 'mult54'
        entry = DefectEntry(defect,bulk_structure,energy_diff,corrections,data,label)
        return entry 
    
    
    def test_formation_energy(self):
        formation_energy = self.entry.formation_energy(self.vbm,self.chempots,fermi_level=0.25)
        self.assert_all_close(formation_energy,4.423740193780213,rtol=100,atol=1e-02)
        
    def test_defect_concentration(self):
        concentration = self.entry.defect_concentration(self.vbm,self.chempots,temperature=300,fermi_level=0.25) 
        self.assert_all_close(concentration,2.4154948187001994e-52,rtol=100,atol=1e-02)
        
        def_conc = self.entry.defect_concentration(self.vbm,self.chempots,temperature=300,fermi_level=5,occupation_function='MB')
        self.assert_all_close(def_conc,1.5117989733276857e+28,rtol=100,atol=1e-02)
        
        def_conc = self.entry.defect_concentration(self.vbm,self.chempots,temperature=300,fermi_level=5,occupation_function='FD') 
        self.assert_all_close(def_conc,4.995386515296238e+22,rtol=100,atol=1e-02)
        
    def test_relaxation_volume(self):
        relaxation_volume = self.entry.relaxation_volume(self.stress_bulk,self.bulk_modulus) # doesn't make physical sense with charged defect
        self.assert_all_close(relaxation_volume,-15.063823947648379)
        
    def test_as_dict_from_dict(self):
        entry_dict_1 = self.entry.as_dict()
        entry_dict_2 = DefectEntry.from_dict(self.entry.as_dict()).as_dict()
        self.assert_object_almost_equal(entry_dict_1,entry_dict_2)