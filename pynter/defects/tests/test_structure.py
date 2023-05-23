#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:01:20 2023

@author: villa
"""

import os 
import os.path as op

from pymatgen.core.composition import Composition

from pynter.data.datasets import Dataset
from pynter.defects.structure import create_def_structure_for_visualization, defect_finder

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/defects/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

ds = Dataset.from_json(get_path('ds_Si-defects.json'))

structure_bulk = ds.select_jobs(names=['Si-bulk_PBE-SCF'])[0].final_structure
structure_vac = ds.select_jobs(names=['Si-test_Vac_Si_q0_PBE-rel_2'])[0].initial_structure
structure_sub = ds.select_jobs(names=['Si-test_Sub_P_on_Si_q1_PBE-rel_2'])[0].initial_structure
structure_int = ds.select_jobs(names=['Si-test_Int_Si(mult54)_q0_PBE-rel_2'])[0].initial_structure

ds_complex = Dataset.from_json(get_path('ds_Si-defects_complexes.json'))

structure_comp = ds_complex.select_jobs(names=['Si-test_Sub_P_on_Si-Sub_B_on_Si_q0_PBE-rel_2'])[0].initial_structure


def test_defect_finder():
    vac = defect_finder(structure_vac,structure_bulk)
    assert vac.defect_type == 'Vacancy'
    assert vac.defect_site_index == 0
    assert vac.defect_specie == 'Si'
    
    sub = defect_finder(structure_sub,structure_bulk)
    assert sub.defect_type == 'Substitution'
    assert sub.defect_site_index == 0
    assert sub.defect_specie == 'P'
    
    inter = defect_finder(structure_int,structure_bulk)
    assert inter.defect_type == 'Interstitial'
    assert inter.defect_site_index == 54
    assert inter.defect_specie == 'Si'
    
    comp = defect_finder(structure_comp,structure_bulk)
    assert comp.defect_type == 'DefectComplex'
    assert comp.defect_names == ['Sub_B_on_Si', 'Sub_P_on_Si']
    assert comp.delta_atoms == {'B': 1, 'Si': -2, 'P': 1}
    
def test_create_def_structure_for_visualization():
    structure_vis = create_def_structure_for_visualization(structure_vac, structure_bulk,sort_to_bulk=True)
    assert structure_vis.composition == Composition('Si53P1')
    assert structure_vis[0].specie.symbol == 'P'
    
    
    
    