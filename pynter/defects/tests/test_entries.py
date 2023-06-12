#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 16:47:57 2023

@author: villa
"""
import os
import os.path as op
from numpy.testing import assert_almost_equal
import unittest

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots
from pynter.defects.defects import Interstitial
from pynter.defects.entries import DefectEntry

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
stress_bulk = [[19.5106052, -0.0, -0.0], [0.0, 19.51060527, 0.0], [0.0, 0.0, 19.51060528]]
bulk_modulus = 144

def get_defect_entry():
    bulk_structure = get_object_from_json(Structure,get_path('Si-bulk_structure_3x3x3_supercell.json'))
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


def test_defect_entry():
    entry = get_defect_entry()
    formation_energy = entry.formation_energy(vbm,chempots,fermi_level=0.25)
    assert_almost_equal(formation_energy,4.423740193780213,3)
    
    concentration = entry.defect_concentration(vbm,chempots,temperature=300,fermi_level=0.25) 
    assert_almost_equal(concentration,2.4154948187001994e-52,3)
    
    def_conc = entry.defect_concentration(vbm,chempots,temperature=300,fermi_level=5,occupation_function='MB')
    assert_almost_equal(def_conc,1.5117989733276857e+28,3)
    
    def_conc = entry.defect_concentration(vbm,chempots,temperature=300,fermi_level=5,occupation_function='FD') 
    assert_almost_equal(def_conc,4.995386515296238e+22,3)
    
    relaxation_volume = entry.relaxation_volume(stress_bulk,bulk_modulus) # doesn't make physical sense with charged defect
    assert_almost_equal(relaxation_volume,-15.063823947648379,3)
    
    entry_dict_1 = entry.as_dict()
    entry_dict_2 = DefectEntry.from_dict(entry.as_dict()).as_dict()
    unittest.TestCase().assertDictEqual(entry_dict_1,entry_dict_2)