#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 16:47:57 2023

@author: villa
"""
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure

from pynter.tools.utils import get_object_from_json
from pynter.defects.chempots.core import Chempots
from pynter.defects.defects import Interstitial, Vacancy
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
        defect = Interstitial(
                            specie = 'Si',
                            defect_site=defect_site,
                            bulk_structure=bulk_structure,
                            charge=-1.0,
                            multiplicity=54)
        
        energy_diff = 5.168792819999965
        corrections = {'kumagai': -0.08825262621975172}
        data = {'stress': 
                [[39.5772376, -0.0, 0.0],
                 [0.0, 39.57723773, 0.0],
                 [0.0, -0.0, 39.57723775]]}
        label = 'mult54'
        entry = DefectEntry(defect,energy_diff,corrections,data,label)
        return entry 
    
    
    def test_formation_energy(self):
        formation_energy = self.entry.formation_energy(self.vbm,self.chempots,fermi_level=0.25)
        self.assert_all_close(formation_energy,4.423740193780213,rtol=100,atol=1e-02)
        
    def test_defect_concentration(self):
        concentration = self.entry.defect_concentration(self.vbm,self.chempots,temperature=300,fermi_level=0.25) 
        self.assert_all_close(concentration,2.4154948187001994e-52,rtol=100,atol=1e-02)
        
        def_conc = self.entry.defect_concentration(self.vbm,self.chempots,temperature=300,fermi_level=5) 
        self.assert_all_close(def_conc,4.995386515296238e+22,rtol=100,atol=1e-02)
        
    def test_relaxation_volume(self):
        relaxation_volume = self.entry.relaxation_volume(self.stress_bulk,self.bulk_modulus) # doesn't make physical sense with charged defect
        self.assert_all_close(relaxation_volume,-15.063823947648379)

    def test_label(self):
        entry = self.entry
        label = 'test'
        entry.set_label(label)
        actual = entry.label
        desired = label
        self.assertEqual(actual, desired)
        
    def test_as_dict_from_dict(self):
        entry_dict_1 = self.entry.as_dict()
        entry_dict_2 = DefectEntry.from_dict(self.entry.as_dict()).as_dict()
        self.assert_object_almost_equal(entry_dict_1,entry_dict_2)


    def test_custom_functions(self):
        entry = DefectEntry(Vacancy('O',charge=2,bulk_volume=800,multiplicity=1),energy_diff=7)
        chempots = {'O': -4.95}
        vbm = 0
        band_gap = 2

        def custom_eform(entry,vbm=None,chemical_potentials=None,fermi_level=0,temperature=300,**kwargs):
                formation_energy = entry.energy_diff + entry.charge*(vbm+fermi_level) 
                if chemical_potentials:
                    chempot_correction = -1 * sum([entry.delta_atoms[el]*chemical_potentials[el] for el in entry.delta_atoms])
                else:
                    chempot_correction = 0
                    
                formation_energy = formation_energy + chempot_correction
                test_quantity = kwargs['test'] if kwargs and 'test' in kwargs.keys() else 0
                return formation_energy - 1/500 *temperature + test_quantity       
        

        def custom_dconc(entry,vbm=0,chemical_potentials=None,temperature=0,fermi_level=0,per_unit_volume=True,eform_kwargs={},**kwargs):
                from pynter.defects.entries import fermi_dirac
                n = entry.defect.site_concentration_in_cm3 if per_unit_volume else entry.multiplicity 
                if per_unit_volume:
                    n = n * (1+1/500*temperature)
                eform = entry.formation_energy(
                                        vbm=vbm,
                                        chemical_potentials=chemical_potentials,
                                        fermi_level=fermi_level,
                                        temperature=temperature,
                                        **eform_kwargs)
                
                test_quantity = kwargs['test_conc'] if kwargs and 'test_conc' in kwargs.keys() else 1
                return  n * fermi_dirac(eform,temperature) * test_quantity
        
        actual = entry.formation_energy(vbm=vbm,chemical_potentials=chempots,fermi_level=1)
        desired = 4.05
        self.assert_all_close(actual, desired)
        actual = entry.formation_energy(vbm=vbm,chemical_potentials=chempots,fermi_level=1,temperature=500)
        self.assert_all_close(actual, desired)

        entry.set_formation_energy_function(custom_eform)
        actual = entry.formation_energy(vbm=vbm,chemical_potentials=chempots,fermi_level=1,temperature=0)
        self.assert_all_close(actual, desired)
        actual = entry.formation_energy(vbm=vbm,chemical_potentials=chempots,fermi_level=1,temperature=500)
        desired = 3.05
        self.assert_all_close(actual, desired)

        actual = entry.formation_energy(vbm=vbm,chemical_potentials=chempots,fermi_level=1,temperature=500,test=2)
        desired = 5.05
        self.assert_all_close(actual, desired)

        actual = entry.formation_energy(vbm=vbm,chemical_potentials=chempots,fermi_level=1,test=2)
        desired = 6.05
        self.assert_all_close(actual, desired)


        actual = entry.defect_concentration(vbm=0,chemical_potentials=chempots,fermi_level=0,temperature=800)
        desired = 1.8255765558702477e+18
        self.assert_all_close(actual, desired, rtol=1e-03)

        entry.set_defect_concentration_function(custom_dconc)
        actual = entry.defect_concentration(vbm=0,chemical_potentials=chempots,fermi_level=0,temperature=800)
        desired = 4.746499045262644e+18
        self.assert_all_close(actual, desired, rtol=1e-03)

        actual = entry.defect_concentration(vbm=0,chemical_potentials=chempots,fermi_level=0,temperature=800,eform_kwargs={'test':1})
        desired = 2383885542239.6626
        self.assert_all_close(actual, desired, rtol=1e-03)

        actual = entry.defect_concentration(vbm=0,chemical_potentials=chempots,fermi_level=0,temperature=800,eform_kwargs={'test':1},test_conc=10)
        desired = desired * 10
        self.assert_all_close(actual, desired, rtol=1e-03)

        actual = entry.defect_concentration(vbm=0,chemical_potentials=chempots,fermi_level=0,temperature=800,test_conc=10)
        desired = 4.746499045262644e+19
        self.assert_all_close(actual, desired, rtol=1e-03)


