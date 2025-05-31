#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 14:33:12 2020

@author: villa
"""

from abc import ABCMeta
from monty.json import MSONable, MontyDecoder
import numpy as np
import warnings

from pymatgen.core.units import kb

from pynter.defects.structure import defect_finder
from pynter.defects.elasticity import Stresses
from pynter.vasp.utils import get_charge_from_computed_entry



class DefectEntry(MSONable,metaclass=ABCMeta):
    
    def __init__(self,defect,energy_diff,corrections,data=None,label=None):
        """
        Contains the data for a defect calculation.
        
        Args:
            defect: Defect object (Vacancy, Interstitial, Substitution, Polaron or DefectComplex)
            energy_diff (float): difference btw energy of defect structure and energy of pure structure
            corrections (dict): Dict of corrections for defect formation energy. All values will be summed and
                                added to the defect formation energy.     
            data : (dict), optional
                Store additional data in dict format.
            label : (str), optional
                Additional label to add to defect specie. Does not influence non equilibrium calculations.
        """
        self._defect = defect
        self._energy_diff = energy_diff
        self._corrections = corrections if corrections else {}
        self._data = data if data else {}
        self._defect.set_label(label)
        
    
    def __repr__(self):
        return "DefectEntry: Name=%s, Charge=%i" %(self.name,self.charge)

    def __str__(self):
        output = [
            "DefectEntry",
            "Defect: %s" %(self.defect.__str__()),
            "Bulk System: %s" %self.bulk_structure.composition,
            "Energy: %.4f" %self.energy_diff,
            "Corrections: %.4f" %sum([v for v in self.corrections.values()]),
            "Charge: %i" %self.charge,
            "Multiplicity: %i" %self.multiplicity,
            "Data: %s" %list(self.data.keys()),
            "Name: %s" %self.name,
            "\n"
            ]
        return "\n".join(output)

    @property
    def bulk_structure(self):
        return self.defect.bulk_structure   

    @property
    def charge(self):
        return self.defect.charge
    
    @property
    def corrections(self):
        return self._corrections
    
    @property
    def data(self):
        """
        Dictionary containing additional data.
        """
        return self._data

    @data.setter
    def data(self,data):
        self._data = data
        return 

    @property
    def defect(self):
        return self._defect  
    
    @property
    def defect_specie(self):
        return self.defect.defect_specie
    
    @property
    def defect_type(self):
        return self.defect.defect_type
    
    @property
    def delta_atoms(self):
        """
        Dictionary with Element as keys and particle difference between defect structure
        and bulk structure as values.
        """
        return self.defect.delta_atoms
    
    @property
    def energy_diff(self):
        return self._energy_diff
    
    @property
    def label(self):
        return self.defect.label
    
    @label.setter
    def label(self,label):
        self.defect.set_label(label)
        return

    @property
    def multiplicity(self):
        return self.defect.multiplicity
    
    @multiplicity.setter
    def multiplicity(self,multiplicity):
        self.defect.set_multiplicity(multiplicity)
        return
            
    @property
    def name(self):
        return self.defect.name
    
    @property
    def symbol(self):
        return self.defect.name.symbol

    @property
    def symbol_charge(self):
        return self.defect.symbol_with_charge

    @property
    def symbol_kroger(self):
        return self.defect.symbol_with_charge_kv
    

    @classmethod
    def from_dict(cls,d):
        """
        Reconstitute a DefectEntry object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of DefectEntry.

        Returns:
            DefectEntry object
        """
        defect = MontyDecoder().process_decoded(d['defect'])
        energy_diff = d['energy_diff']
        corrections = d['corrections']
        data = d['data']
        label = d['label']
        return cls(defect=defect,energy_diff=energy_diff,corrections=corrections,data=data,label=label)

    @staticmethod
    def from_computed_entries(computed_entry_defect,computed_entry_bulk,corrections,
                              multiplicity=1,data=None,label=None,tol=1e-03):
        """
        Generate DefectEntry object from Pymatgen's ComputedStructureEntry objects.

        Parameters
        ----------
        computed_entry_defect : (VaspJob)
            ComputedStructureEntry of the defect calculation.
        job_bulk : (VaspJob)
            ComputedStructureEntry of the bulk calculation.
        corrections : (dict)
            Dict of corrections for defect formation energy. All values will be summed and
            added to the defect formation energy.
        multiplicity : (int), optional
            Multiplicity of defect within the supercell. 
            If set to None is attempted to be determined automatically with Pymatgen. The default is 1.
        data : (dict), optional
            Store additional data in dict format.
        label : (str), optional
            Additional label to add to defect specie. Does not influence non equilibrium calculations.
        tol : (float)
            Tolerance for defect_finder function. The default is 1e-03.

        Returns
        -------
        DefectEntry
        """ 
        entry_df, entry_bulk = computed_entry_defect,computed_entry_bulk
        charge = get_charge_from_computed_entry(entry_df)
        energy_diff = entry_df.energy - entry_bulk.energy
        
        return DefectEntry.from_structures(entry_df.structure, entry_bulk.structure, energy_diff,
                                           corrections,charge,multiplicity,data,label,tol=tol)
        

    @staticmethod
    def from_jobs(job_defect, job_bulk, corrections, defect_structure=None,
                  multiplicity=1,data=None,label=None,tol=1e-03):
        """
        Generate DefectEntry object from VaspJob objects.

        Parameters
        ----------
        job_defect : (VaspJob)
            Defect calculation.
        job_bulk : (VaspJob)
            Bulk calculation.
        corrections : (dict)
            Dict of corrections for defect formation energy. All values will be summed and
            added to the defect formation energy.
        defect_structure : (Structure)
            Structure of the defect. If None the intial structure of job_defect is taken. The default is None. 
        multiplicity : (int), optional
            Multiplicity of defect within the supercell. 
            If set to None is attempted to be determined automatically with Pymatgen. The default is 1.
        data : (dict), optional
            Store additional data in dict format.
        label : (str), optional
            Additional label to add to defect specie. Does not influence non equilibrium calculations.
        tol : (float)
            Tolerance for defect_finder function. The default is 1e-03.

        Returns
        -------
        DefectEntry
        """ 
        defect_structure = defect_structure if defect_structure else job_defect.initial_structure
        bulk_structure = job_bulk.final_structure
        energy_diff = job_defect.final_energy - job_bulk.final_energy
        charge = job_defect.charge
        
        return DefectEntry.from_structures(defect_structure, bulk_structure, energy_diff, corrections,
                                                  charge,multiplicity,data,label,tol=tol)


    @staticmethod
    def from_structures(defect_structure,bulk_structure,energy_diff,corrections,charge=0,
                        multiplicity=1,data=None,label=None,tol=1e-03):
        """
        Generate DefectEntry object from Structure objects.

        Parameters
        ----------
        defect_structure : (Structure)
            Defect structure.
        bulk_structure : (Structure)
            Bulk structure.
        energy_diff (float): 
            Difference btw energy of defect structure and energy of pure structure
        corrections (dict): 
            Dict of corrections for defect formation energy. All values will be summed and
            added to the defect formation energy.  
        charge : (int), optional
            Charge of the defect system. The default is 0.
        multiplicity : (int), optional
            multiplicity of defect within the supercell. 
            If set to None is attempted to be determined automatically with Pymatgen. The default is 1.
        data : (dict), optional
            Store additional data in dict format.
        label : (str), optional
            Additional label to add to defect specie. Does not influence non equilibrium calculations.
        tol : (float)
            Tolerance for defect_finder function. The default is 1e-03.

        Returns
        -------
        DefectEntry
        """
        defect = defect_finder(defect_structure, bulk_structure,tol=tol)
        defect.set_charge(charge)
        if multiplicity:
            defect.set_multiplicity(multiplicity)
        else:
            try:
                new_multiplicity = defect.get_multiplicity()
                defect.set_multiplicity(new_multiplicity)
            except NotImplementedError:
                warnings.warn(f'get_multiplicity not implemented for {defect.defect_type}, setting multiplicity to 1')
                defect.set_multiplicity(1)
        
        return DefectEntry(defect, energy_diff, corrections,data,label)


    def defect_concentration(self, vbm, chemical_potentials, temperature=300, fermi_level=0.0, 
                             per_unit_volume=True,occupation_function='MB'):
        """
        Compute the defect concentration for a temperature and Fermi level.
        Args:
            temperature:
                the temperature in K
            fermi_level:
                the fermi level in eV (with respect to the VBM)
        Returns:
            defects concentration in cm^-3
        """
        n = self.defect.site_concentration_in_cm3 if per_unit_volume else self.multiplicity 
        eform = self.formation_energy(vbm, chemical_potentials, fermi_level=fermi_level)
        
        if occupation_function=='FD':
            conc = n * fermi_dirac(eform,temperature)
        elif occupation_function=='MB':
            conc = n * maxwell_boltzmann(eform,temperature)
        else:
            raise ValueError('Invalid occupation function. Options are: "FD" for Fermi-Dirac and "MB" for Maxwell-Boltzmann.')
        return conc


    def formation_energy(self,vbm,chemical_potentials,fermi_level=0):
        """
        Compute the formation energy for a defect taking into account a given chemical potential and fermi_level
        Args:
            vbm(float): Valence band maximum of pure structure
            chemical_potentials (dict): Dictionary of elemental chemical potential values.
                Keys are Element objects within the defect structure's composition.
                Values are float numbers equal to the atomic chemical potential for that element.
            fermi_level (float):  Value corresponding to the electron chemical potential.
            """
            
        formation_energy = (self.energy_diff + self.charge*(vbm+fermi_level) + 
                       sum([ self.corrections[correction_type]  for correction_type in self.corrections ]) 
                        ) 
        
        if chemical_potentials:
            chempot_correction = -1 * sum([self.delta_atoms[el]*chemical_potentials[el] for el in self.delta_atoms])
        else:
            chempot_correction = 0
            
        formation_energy = formation_energy + chempot_correction
        
        return formation_energy
    
    
    def relaxation_volume(self,stress_bulk,bulk_modulus,add_corrections=True): #still to decide weather to keep this method
        """
        Calculate relaxation volume from stresses. Stresses data needs to be in numpy.array format and present 
        in the "data" dictionary with realtive "stress" key. Duplicate of function that can be found in Stresses
        class in elasticity module, added here for convenience.

        Parameters
        ----------
        stress_bulk : (np.array)
            Stresses of bulk calculation.
        bulk_modulus : (float)
            Bulk modulus in GPa.
        add_corrections : (bool)
            Add correction terms from "elastic_corrections" dict (if key is present in dict).

        Returns
        -------
        rel_volume : (float)
            Relaxation volume in A°^3.
        """
        es = Stresses(stress_bulk)
        return es.get_relaxation_volume(self, bulk_modulus)



def fermi_dirac(E,T):
    """
    Returns the defect occupation as a function of the formation energy,
    using the Fermi-Dirac distribution with chemical potential equal to 0. 
    Args:
        E (float): energy in eV
        T (float): the temperature in kelvin
    """
    return 1. / (1. + np.exp(E/(kb*T)) )


def maxwell_boltzmann(E,T):
    """
    Returns the defect occupation as a function of the formation energy,
    using the exponential dependence of the Maxwell-Boltzmann distribution. 
    This is the more common approach, which is the approximation of the FD-like 
    distribution for N_sites >> N_defects.
    Args:
        E (float): energy in eV
        T (float): the temperature in kelvin
    """
    return np.exp(-1.0*E /(kb*T)) 
