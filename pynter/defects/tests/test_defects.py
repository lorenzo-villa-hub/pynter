#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 11:59:40 2023

@author: villa
"""
import numpy as np
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.composition import Composition


from pynter.defects.defects import Vacancy, Interstitial, Substitution, Polaron, DefectComplex

from pynter.testing.core import PynterTest
from pynter.testing.defects import DefectTest


bulk_structure = PynterTest().structure.copy()
bulk_structure.make_supercell(3)


class TestDefect(PynterTest):

    def test_interstitial(self):
        structure = bulk_structure.copy()
        frac_coords = np.array([0.66666667, 0.5, 0.5])
        interstitial_site = PeriodicSite('Si', frac_coords, structure.lattice)
        inter = Interstitial(
                           specie='Si',
                           defect_site=interstitial_site,
                           bulk_structure=structure,
                           charge=0,
                           multiplicity=108,
                           label='test') 
        
        assert inter.specie == 'Si'
        assert inter.type == 'Interstitial'
        assert inter.defect_site_index == 54
        assert inter.defect_composition == Composition('Si55')
        assert inter.defect_structure.composition == inter.defect_composition
        assert inter.name == 'Int_Si(test)'
        assert inter.charge == 0
        assert inter.delta_atoms == {'Si': 1}
        assert inter.symbol == '$Si_i$(test)'
        assert inter.symbol_with_charge == '$Si_i$(test)$^{\\;0}$'
        assert inter.symbol_with_charge_kv == '$Si_i$(test)$^{x}$'
    
    
    def test_polaron(self):
        structure = bulk_structure.copy()
        site = structure[0]
        pol = Polaron(
                    specie='Si',
                    defect_site=site,
                    bulk_structure=structure,
                    charge=1,
                    multiplicity=1,
                    label='test') 

        assert pol.specie == 'Si'
        assert pol.type == 'Polaron'
        assert pol.defect_site_index == 0
        assert pol.defect_composition == Composition('Si54')
        assert pol.defect_structure.composition == pol.defect_composition
        assert pol.charge == 1
        assert pol.delta_atoms == {}
        assert pol.get_multiplicity() == 54
        pol.set_multiplicity(54)
        assert pol.multiplicity == 54
        assert pol.name == 'Pol_Si(test)'
        assert pol.symbol == '$Si_{Si}$(test)'
        assert pol.symbol_with_charge =='$Si_{Si}$(test)$^{+1}$'
        assert pol.symbol_with_charge_kv == '$Si_{Si}$(test)$^{°}$'
    
    
    def test_substitution(self):
        structure = bulk_structure.copy()
        site = structure[0]
        defect_site = PeriodicSite('P',site.frac_coords,site.lattice)
        sub = Substitution(
                    specie='P',
                    bulk_specie='Si',
                    defect_site=site,
                    bulk_structure=structure,
                    charge=1,
                    multiplicity=1,
                    label='test') 
 
        assert sub.specie == 'P'
        assert sub.type == 'Substitution'
        assert sub.defect_site_index == 0
        assert sub.defect_composition == Composition('Si53P1')
        assert sub.defect_structure.composition == sub.defect_composition
        assert sub.name == 'Sub_P_on_Si(test)'
        assert sub.charge == 1
        assert sub.delta_atoms == {'P': 1, 'Si': -1}
        assert sub.site_in_bulk == site
        assert sub.get_multiplicity() == 54
        sub.set_multiplicity(54)
        assert sub.multiplicity == 54
        assert sub.symbol == '$P_{Si}$(test)'
        assert sub.symbol_with_charge == '$P_{Si}$(test)$^{+1}$'
        assert sub.symbol_with_charge_kv == '$P_{Si}$(test)$^{°}$'
        
    
    def test_vacancy(self):
        structure = bulk_structure.copy()
        site = structure[0]
        vac = Vacancy(
                           specie='Si',
                           defect_site=site,
                           bulk_structure=structure,
                           charge=0,
                           multiplicity=1,
                           label='test') 
        vac = Vacancy(site,structure,charge=0,multiplicity=1,label='test')
        assert vac.specie == 'Si'
        assert vac.type == 'Vacancy'
        assert vac.defect_site_index == 0
        assert vac.defect_composition == Composition('Si53')
        assert vac.defect_structure.composition == vac.defect_composition
        assert vac.charge == 0
        assert vac.delta_atoms == {'Si':-1}
        assert vac.get_multiplicity() == 54
        vac.set_multiplicity(54)
        assert vac.multiplicity == 54
        assert vac.name == 'Vac_Si(test)'
        assert vac.symbol == '$V_{Si}$(test)'
        assert vac.symbol_with_charge =='$V_{Si}$(test)$^{\\;0}$'
        assert vac.symbol_with_charge_kv == '$V_{Si}$(test)$^{x}$'
        
        
    def test_defect_complex(self):
        structure = bulk_structure.copy()
        vac_site = structure[0]
        vac = Vacancy(
                           specie='Si',
                           defect_site=vac_site,
                           bulk_structure=structure,
                           charge=0,
                           multiplicity=54,
                           label='vac') 

        sub_site = structure[1]
        defect_site = PeriodicSite('P',sub_site.frac_coords,sub_site.lattice)
        sub = Substitution(
                           specie='P',
                           bulk_specie='Si',
                           defect_site=sub_site,
                           bulk_structure=structure,
                           charge=1,
                           multiplicity=54,
                           label='sub') 

        defects = [vac,sub]
        comp = DefectComplex(
                            defects=defects,
                            bulk_structure=structure,
                            multiplicity=54*4,
                            charge=1,
                            label='test')
        assert comp.type == 'DefectComplex'
        assert comp.name == 'Vac_Si(vac)-Sub_P_on_Si(sub)(test)'
        assert comp.defect_composition == Composition('Si52P1')
        assert comp.charge == 1
        assert comp.delta_atoms == {'Si': -2, 'P': 1}
        assert comp.multiplicity == 216
        assert comp.symbol == '$V_{Si}$-$P_{Si}$(test)'
        assert comp.symbol_with_charge == '$V_{Si}$-$P_{Si}$(test)$^{+1}$'
        assert comp.symbol_with_charge_kv == '$V_{Si}$-$P_{Si}$(test)$^{°}$'
    
    
    
    
# unit_structure = PynterTest().structure.copy()

# def test_create_vacancies():
#     vac = create_vacancies(unit_structure,['Si'],supercell_size=3)[0]
#     structure = unit_structure.copy()
#     structure.make_supercell(3)
#     frac_coords = structure[0].frac_coords
#     defect_site = PeriodicSite('Si', frac_coords, structure.lattice)
#     vac_test = Vacancy(defect_site, structure)
#     DefectTest().assert_Defect_equal(vac,vac_test)


# def test_create_substitutions():    
#     sub = create_substitutions(unit_structure,elements_to_replace={'Si':'P'},supercell_size=3)[0]
#     structure = unit_structure.copy()
#     structure.make_supercell(3)
#     frac_coords = structure[0].frac_coords
#     defect_site = PeriodicSite('P', frac_coords, structure.lattice)
#     sub_test = Substitution(defect_site, structure)
#     DefectTest().assert_Defect_equal(sub,sub_test)
    

# commented out because of long run time, uncomment to include in tests
# def test_create_interstitials():
#     inter = create_interstitials(unit_structure,['Si'],supercell_size=3)[0]
#     structure = bulk_structure.copy()
#     frac_coords = np.array([0.66666667, 0.5, 0.5])
#     interstitial_site = PeriodicSite('Si', frac_coords, structure.lattice)
#     multiplicity = 108
#     inter_test = Interstitial(interstitial_site,structure,multiplicity=108,label='mult108')   
#     CompareDefects().compare(inter,inter_test)
    




















