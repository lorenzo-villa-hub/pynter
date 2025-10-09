#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:01:20 2023

@author: villa
"""
import os.path as op

from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.io.vasp.outputs import Vasprun

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
        assert vac.type == 'Vacancy'
        assert vac.defect_site_index == 0
        assert vac.specie == 'Si'
        
    def test_finder_substitution(self):
        sub = defect_finder(self.structure_sub,self.structure_bulk)
        assert sub.type == 'Substitution'
        assert sub.defect_site_index == 0
        assert sub.specie == 'P'
        
    def test_finder_interstitial(self):
        inter = defect_finder(self.structure_int,self.structure_bulk)
        assert inter.type == 'Interstitial'
        assert inter.defect_site_index == 54
        assert inter.specie == 'Si'
        
    def test_finder_complex(self):
        comp = defect_finder(self.structure_comp,self.structure_bulk)
        assert comp.type == 'DefectComplex'
        assert comp.defect_names == ['Sub_B_on_Si', 'Sub_P_on_Si']
        assert comp.delta_atoms == {'B': 1, 'Si': -2, 'P': 1}

    def test_finder_max_defects(self):
        structure_problematic = Vasprun(
            self.get_testfile_path('SiO2-defects/Defects/vacancies/Vac_O/q2/2-PBE-OPT/vasprun.xml'),
            parse_dos=False,parse_eigen=False,parse_potcar_file=False).final_structure
        structure_bulk = Vasprun(
            self.get_testfile_path('SiO2-defects/Bulk-2x2x2-supercell/vasprun.xml'),
            parse_dos=False,parse_eigen=False,parse_potcar_file=False).final_structure
        defect = defect_finder(
                            structure_defect=structure_problematic,
                            structure_bulk=structure_bulk,
                            tol=1e-02,
                            max_number_of_defects=1)
        assert defect.type == 'Vacancy'
        assert defect.specie == 'O'
        assert defect.defect_site_index == 24

        
    def test_create_def_structure_for_visualization(self):
        structure_vis = create_def_structure_for_visualization(self.structure_vac, self.structure_bulk,sort_to_bulk=True)
        assert structure_vis.composition == Composition('Si53P1')
        assert structure_vis[0].specie.symbol == 'P'
    
    
    
    