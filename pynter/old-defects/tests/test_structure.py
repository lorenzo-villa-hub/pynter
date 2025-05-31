#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:01:20 2023

@author: villa
"""
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition

from pynter.defects.structure import create_def_structure_for_visualization, defect_finder
from pynter.tools.utils import get_object_from_json

from pynter.testing.core import PynterTest


class TestDefectFinder(PynterTest):

    def setUp(self):
        self.structure_bulk = get_object_from_json(Structure,self.get_testfile_path('structure_bulk.json'))
        self.structure_vac = get_object_from_json(Structure,self.get_testfile_path('structure_vac.json'))
        self.structure_sub = get_object_from_json(Structure,self.get_testfile_path('structure_sub.json'))
        self.structure_int = get_object_from_json(Structure,self.get_testfile_path('structure_int.json'))
        self.structure_comp = get_object_from_json(Structure,self.get_testfile_path('structure_comp.json'))

    def test_finder_vacancy(self):
        vac = defect_finder(self.structure_vac,self.structure_bulk)
        assert vac.defect_type == 'Vacancy'
        assert vac.defect_site_index == 0
        assert vac.defect_specie == 'Si'
        
    def test_finder_substitution(self):
        sub = defect_finder(self.structure_sub,self.structure_bulk)
        assert sub.defect_type == 'Substitution'
        assert sub.defect_site_index == 0
        assert sub.defect_specie == 'P'
        
    def test_finder_interstitial(self):
        inter = defect_finder(self.structure_int,self.structure_bulk)
        assert inter.defect_type == 'Interstitial'
        assert inter.defect_site_index == 54
        assert inter.defect_specie == 'Si'
        
    def test_finder_complex(self):
        comp = defect_finder(self.structure_comp,self.structure_bulk)
        assert comp.defect_type == 'DefectComplex'
        assert comp.defect_names == ['Sub_B_on_Si', 'Sub_P_on_Si']
        assert comp.delta_atoms == {'B': 1, 'Si': -2, 'P': 1}
        
    def test_create_def_structure_for_visualization(self):
        structure_vis = create_def_structure_for_visualization(self.structure_vac, self.structure_bulk,sort_to_bulk=True)
        assert structure_vis.composition == Composition('Si53P1')
        assert structure_vis[0].specie.symbol == 'P'
    
    
    
    