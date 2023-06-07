#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:02:03 2023

@author: villa
"""
import os
import os.path as op
from numpy.testing import assert_almost_equal

from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots, Reservoirs, PDHandler


homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/phase_diagram/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

mu_refs = {'Na': -1.31, 'Nb': -10.1, 'O': -4.95}
chempots_dict = {'Na':-2.26,'Nb':-6.15,'O':-1.92}
mu = Chempots(chempots_dict)
mu_abs = Chempots({'Na':-3.57,'Nb':-16.25,'O':-6.87})

def test_chempots():
    assert mu == Chempots.from_dict(mu.as_dict())
    mu_elements = {Element('Na'):-2.26,Element('Nb'):-6.15,Element('O'):-1.92}
    assert mu == Chempots.from_pmg_elements(mu_elements)
    assert mu.to_pmg_elements() == mu_elements
    assert mu.get_absolute(mu_refs) == mu_abs
    assert mu_abs.get_referenced(mu_refs) == mu
    
def test_reservoirs():
    res_dict = {'X':mu}
    pd = get_object_from_json(PhaseDiagram,get_path('PD_Na-Nb-O.json'))
    res = Reservoirs(res_dict,pd,are_chempots_delta=True)
    assert res.mu_refs == Chempots(mu_refs)
    assert res == Reservoirs.from_dict(res.as_dict())
    assert res == Reservoirs.from_json(get_path('reservoirs_NaNbO3.json'))
    res.set_to_absolute()
    assert res == Reservoirs({'X':mu_abs})
    res.set_to_referenced()
    assert res == Reservoirs({'X':mu})
    res_filtered = Reservoirs({'X':Chempots({'Na':-2.26,'O':-1.92})})
    assert res.filter_reservoirs(elements=['Na','O']) == res_filtered
    
def test_pdhandler():
    comp = Composition('NaNbO3')
    pd = get_object_from_json(PhaseDiagram,get_path('PD_Na-Nb-O.json'))
    pdh = PDHandler(pd)
    assert_almost_equal(pdh.calculate_single_chempot(comp, Chempots({'O':-1.92,'Na':-2.26})),-6.26,decimal=2)
    boundary_chempots = {
     'NaNbO3-Na2Nb3O6-Nb12O29': {'O': -8.95, 'Nb': -10.86, 'Na': -2.84},
     'NaNbO3-Na3NbO4-O2': {'O': -4.95, 'Nb': -21.43, 'Na': -4.26},
     'NaNbO3-Na2Nb3O6-NaNbO2': {'O': -9.08, 'Nb': -10.86, 'Na': -2.44},
     'NaNbO3-Na3NbO4-NaNbO2': {'O': -9.08, 'Nb': -11.1, 'Na': -2.2},
     'NaNb3O8-NaNbO3-Nb12O29': {'O': -8.14, 'Nb': -12.81, 'Na': -3.31},
     'NaNb3O8-NaNbO3-O2': {'O': -4.95, 'Nb': -20.79, 'Na': -4.91}
     }
    assert pdh.get_all_boundaries_chempots(comp) == boundary_chempots
    target_entries=[]
    for e in pd.all_entries:
        if e.composition.reduced_composition == comp:
            target_entries.append(e)
    assert pdh.get_entries_from_comp(comp) == target_entries
    assert pdh.get_stable_entry_from_comp(comp) == pd.entries[94]
    assert_almost_equal (pdh.get_formation_energy_from_stable_comp(comp), -14.28,decimal=2)
    phase_boundary_chempots = {
     'NaNb3O8-NaNbO3': {'Na': -2.64, 'Nb': -5.89, 'O': -1.92},
     'NaNbO3-Na3NbO4': {'Na': -1.99, 'Nb': -6.53, 'O': -1.92}
     }
    assert pdh.get_phase_boundaries_chempots(comp,{'O':-1.92}) == phase_boundary_chempots
    
    
    
    
    
    
    
    
    
    