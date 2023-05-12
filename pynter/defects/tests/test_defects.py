#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 11:59:40 2023

@author: villa
"""
import unittest
import numpy as np
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.composition import Composition

from pynter.tests.__init__ import get_structure_Si
from pynter.defects.defects import *
from pynter.defects.tests.compare import CompareDefects


bulk_structure = get_structure_Si().copy()
bulk_structure.make_supercell(3)

    

def test_interstitial():
    structure = bulk_structure.copy()
    frac_coords = np.array([0.66666667, 0.5, 0.5])
    interstitial_site = PeriodicSite('Si', frac_coords, structure.lattice)
    multiplicity = 108
    inter = Interstitial(interstitial_site,structure,charge=0,multiplicity=108,label='test')
    assert inter.defect_specie == 'Si'
    assert inter.defect_type == 'Interstitial'
    assert inter.defect_site_index == 54
    assert inter.defect_composition == Composition('Si55')
    assert inter.name == DefectName('Interstitial','Si',label='test')
    assert inter.charge == 0
    assert inter.delta_atoms == {'Si': 1}
    assert inter.symbol == '$Si_{i}$(test)'
    assert inter.symbol_with_charge == '$Si_{i}$(test)$^{\\;0}$'
    assert inter.symbol_with_charge_kv == '$Si_{i}$(test)$^{x}$'


def test_polaron():
    structure = bulk_structure.copy()
    site = structure[0]
    pol = Polaron(site,structure,charge=1,multiplicity=1,label='test')    
    assert pol.defect_specie == 'Si'
    assert pol.defect_type == 'Polaron'
    assert pol.defect_site_index == 0
    assert pol.defect_composition == Composition('Si54')
    assert pol.charge == 1
    assert pol.delta_atoms == {}
    assert pol.get_multiplicity() == 54
    pol.set_multiplicity(54)
    assert pol.multiplicity == 54
    assert pol.name == DefectName('Polaron','Si',label='test')
    assert pol.symbol == '$Si_{Si}$(test)'
    assert pol.symbol_with_charge =='$Si_{Si}$(test)$^{+1}$'
    assert pol.symbol_with_charge_kv == '$Si_{Si}$(test)$^{°}$'


def test_substitution():
    structure = bulk_structure.copy()
    site = structure[0]
    defect_site = PeriodicSite('P',site.frac_coords,site.lattice)
    sub = Substitution(defect_site,structure,charge=1,multiplicity=1,label='test')
    assert sub.defect_specie == 'P'
    assert sub.defect_type == 'Substitution'
    assert sub.defect_site_index == 0
    assert sub.defect_composition == Composition('Si53P1')
    assert sub.name == DefectName('Substitution','P',bulk_specie='Si',label='test')
    assert sub.charge == 1
    assert sub.delta_atoms == {'P': 1, 'Si': -1}
    assert sub.site_in_bulk == site
    assert sub.get_multiplicity() == 54
    sub.set_multiplicity(54)
    assert sub.multiplicity == 54
    assert sub.symbol == '$P_{Si}$(test)'
    assert sub.symbol_with_charge == '$P_{Si}$(test)$^{+1}$'
    assert sub.symbol_with_charge_kv == '$P_{Si}$(test)$^{°}$'
    

def test_vacancy():
    structure = bulk_structure.copy()
    site = structure[0]
    vac = Vacancy(site,structure,charge=0,multiplicity=1,label='test')
    assert vac.defect_specie == 'Si'
    assert vac.defect_type == 'Vacancy'
    assert vac.defect_site_index == 0
    assert vac.defect_composition == Composition('Si53')
    assert vac.charge == 0
    assert vac.delta_atoms == {'Si':-1}
    assert vac.get_multiplicity() == 54
    vac.set_multiplicity(54)
    assert vac.multiplicity == 54
    assert vac.name == DefectName('Vacancy','Si',label='test')
    assert vac.symbol == '$V_{Si}$(test)'
    assert vac.symbol_with_charge =='$V_{Si}$(test)$^{\\;0}$'
    assert vac.symbol_with_charge_kv == '$V_{Si}$(test)$^{x}$'
    
    
def test_defect_complex():
    structure = bulk_structure.copy()
    vac_site = structure[0]
    vac = Vacancy(vac_site,structure,charge=0,multiplicity=54,label='test')
    sub_site = structure[0]
    defect_site = PeriodicSite('P',sub_site.frac_coords,sub_site.lattice)
    sub = Substitution(defect_site,structure,charge=1,multiplicity=54,label='test')
    defects = [vac,sub]
    comp = DefectComplex(defects, structure,multiplicity=54*4,charge=1,label='test')
    assert comp.defect_type == 'DefectComplex'
    assert comp.name == DefectComplexName([DefectName('Vacancy','Si',label='test'),
                                 DefectName('Substitution','P','Si',label='test')],label='test')
    assert comp.defect_names == [DefectName('Vacancy','Si',label='test'),
                                 DefectName('Substitution','P','Si',label='test')]
    assert comp.charge == 1
    assert comp.delta_atoms == {'Si': -2, 'P': 1}
    assert comp.multiplicity == 216
    assert comp.symbol == '$V_{Si}$-$P_{Si}$(test)'
    assert comp.symbol_with_charge == '$V_{Si}$-$P_{Si}$(test)$^{+1}$'
    assert comp.symbol_with_charge_kv == '$V_{Si}$-$P_{Si}$(test)$^{°}$'
    
    
def test_defect_name():
    name = DefectName('Vacancy','Si',label='test')
    assert name.name == 'Vac_Si'
    assert name.label == 'test'
    assert name.fullname == 'Vac_Si(test)'
    assert name.symbol == '$V_{Si}$(test)'
    assert name.dspecie == 'Si'
    assert name.dtype == 'Vacancy'
    assert name == 'Vac_Si(test)'
    assert name == DefectName('Vacancy','Si',label='test')
    assert {name:1} == {'Vac_Si(test)':1}
    assert name == DefectName.from_string('Vac_Si(test)')
    
    name = DefectName('Substitution','P','Si',label='test')
    assert name.name == 'Sub_P_on_Si'
    assert name.label == 'test'
    assert name.fullname == 'Sub_P_on_Si(test)'
    assert name.symbol == '$P_{Si}$(test)'
    assert name.dspecie == 'P'
    assert name.dtype == 'Substitution'
    assert name == 'Sub_P_on_Si(test)'
    assert name == DefectName('Substitution','P','Si',label='test')
    assert {name:1} == {'Sub_P_on_Si(test)':1}
    assert name == DefectName.from_string('Sub_P_on_Si(test)')
    
    name = DefectName('Interstitial','Si',label='test')
    assert name.name == 'Int_Si'
    assert name.label == 'test'
    assert name.fullname == 'Int_Si(test)'
    assert name.symbol == '$Si_{i}$(test)'
    assert name.dspecie == 'Si'
    assert name.dtype == 'Interstitial'
    assert name == 'Int_Si(test)'
    assert name == DefectName('Interstitial','Si',label='test')
    assert {name:1} == {'Int_Si(test)':1}
    assert name == DefectName.from_string('Int_Si(test)')
    
    name = DefectName('Polaron','Si',label='test')
    assert name.name == 'Pol_Si'
    assert name.label == 'test'
    assert name.fullname == 'Pol_Si(test)'
    assert name.symbol == '$Si_{Si}$(test)'
    assert name.dspecie == 'Si'
    assert name.dtype == 'Polaron'
    assert name == 'Pol_Si(test)'
    assert name == DefectName('Polaron','Si',label='test')
    assert {name:1} == {'Pol_Si(test)':1}
    assert name == DefectName.from_string('Pol_Si(test)')


def test_defect_complex_name():
    vac_name = DefectName('Vacancy','Si',label='test')
    sub_name = DefectName('Substitution','P','Si',label='test')
    name = DefectComplexName([vac_name,sub_name],label='test')
    assert name.name == 'Vac_Si-Sub_P_on_Si'
    assert name.label == 'test'
    assert name.fullname =='Vac_Si-Sub_P_on_Si(test)'
    assert name.symbol == '$V_{Si}$-$P_{Si}$(test)'
    assert name.dspecies == ['Si', 'P']
    assert name.dtype == 'DefectComplex'
    assert name == 'Vac_Si-Sub_P_on_Si(test)'
    assert name == DefectComplexName.from_string('Vac_Si-Sub_P_on_Si(test)')
    
    
unit_structure = get_structure_Si().copy()

def test_create_vacancies():
    vac = create_vacancies(unit_structure,['Si'],supercell_size=3)[0]
    structure = unit_structure.copy()
    structure.make_supercell(3)
    frac_coords = structure[0].frac_coords
    defect_site = PeriodicSite('Si', frac_coords, structure.lattice)
    vac_test = Vacancy(defect_site, structure)
    CompareDefects().compare(vac,vac_test)


def test_create_substitutions():    
    sub = create_substitutions(unit_structure,elements_to_replace={'Si':'P'},supercell_size=3)[0]
    structure = unit_structure.copy()
    structure.make_supercell(3)
    frac_coords = structure[0].frac_coords
    defect_site = PeriodicSite('P', frac_coords, structure.lattice)
    sub_test = Substitution(defect_site, structure)
    CompareDefects().compare(sub,sub_test)
    

def test_create_interstitials():
    inter = create_interstitials(unit_structure,['Si'],supercell_size=3)[0]
    structure = bulk_structure.copy()
    frac_coords = np.array([0.66666667, 0.5, 0.5])
    interstitial_site = PeriodicSite('Si', frac_coords, structure.lattice)
    multiplicity = 108
    inter_test = Interstitial(interstitial_site,structure,multiplicity=108,label='mult108')   
    CompareDefects().compare(inter,inter_test)
    




















