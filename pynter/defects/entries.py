#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 14:33:12 2020

@author: villa
"""


import numpy as np
import importlib
from pymatgen.core.units import kb
from pymatgen.core.structure import Structure
from pymatgen.analysis.defects.core import *
from pynter.defects.utils import *
from monty.json import MontyDecoder


def get_defect_entry_from_jobs(job_defect,job_bulk,corrections,defect_structure=None,tol=1e-03,multiplicity=None,label=None):
    """
    Get defect entry from VaspJob objects of defect and bulk calculation.
    The defect_finder method is used to determine if the defect site is single or 
    a list of defect sites (defect complex). In this way the SingleDefectEntry or the 
    DefectComplexEntry is initialized.

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
        Multiplicity of defect within the supercell. The default is None.
        If not provided, for single defects it is determined by pymatgen with symmetry analysis,
        for a defect complex is set to 1.
    tol : (float)
        Tolerance for defect_finder function. The default is 1e-03.
    label : (str), optional
        Additional label to add to defect specie. Does not influence non equilibrium calculations.

    Returns
    -------
    SingleDefectEntry or DefectComplexEntry

    """
    defect_structure = defect_structure if defect_structure else job_defect.initial_structure
    bulk_structure = job_bulk.final_structure
    defects = defect_finder(defect_structure,bulk_structure,tol)
    if isinstance(defects,list):
        entry = DefectComplexEntry.from_jobs(job_defect,job_bulk,corrections,defect_structure,multiplicity,label)
    else:
        entry = SingleDefectEntry.from_jobs(job_defect,job_bulk,corrections,defect_structure,multiplicity,label)
        
    return entry


def get_formatted_legend(fullname):
    # handling label case
    if '(' in fullname:
        fullname = fullname.split('(')
        name = fullname[0]
        entry_label = '('+fullname[1]
    else:
        name = fullname
        entry_label = ''
    # single defect    
    if '-' not in [c for c in name]:        
        flds = name.split('_')
        if 'Vac' == flds[0]:
            base = '$V'
            sub_str = '_{' + flds[1] + '}$'
        elif 'Sub' == flds[0]:
            flds = name.split('_')
            base = '$' + flds[1]
            sub_str = '_{' + flds[3] + '}$'
        elif 'Int' == flds[0]:
            base = '$' + flds[1]
            sub_str = '_{int}$'
        else:
            base = name
            sub_str = ''

        return  base + sub_str + entry_label
    # defect complex
    else:
        label = ''
        names = name.split('-')
        for name in names:
            flds = name.split('_')
            if '-' not in flds:
                if 'Vac' == flds[0]:
                    base = '$V'
                    sub_str = '_{' + flds[1] + '}$'
                elif 'Sub' == flds[0]:
                    flds = name.split('_')
                    base = '$' + flds[1]
                    sub_str = '_{' + flds[3] + '}$'
                elif 'Int' == flds[0]:
                    base = '$' + flds[1]
                    sub_str = '_{int}$'
                else:
                    base = name
                    sub_str = ''
        
                if names.index(name) != (len(names) - 1):
                    label += base + sub_str + '-'
                else:
                    label += base + sub_str

        return label + entry_label


class SingleDefectEntry:
    
    def __init__(self,defect,bulk_structure,energy_diff,corrections,label=None):
        """
        Initializes the data of a single defect. Inspired from the DefectEntry class in pymatgen (pmg).

        Args:
            defect: Pymatgen defect object (Vacancy, Interstitial or Substitution)
            bulk_structure: Pymatgen Structure without any defects
            energy_diff (float): difference btw energy of defect structure and energy of pure structure
            corrections (dict): Dict of corrections for defect formation energy. All values will be summed and
                                added to the defect formation energy.     
            label : (str), optional
                Additional label to add to defect specie. Does not influence non equilibrium calculations.
        """        
        self._defect = defect 
        self._bulk_structure = bulk_structure
        self._energy_diff = energy_diff
        self._corrections = corrections if corrections else {}
        self._label = label        
        

    @property
    def defect(self):
        return self._defect
        
    @property
    def bulk_structure(self):
        return self._bulk_structure
    
    @property
    def energy_diff(self):
        return self._energy_diff
    
    @property
    def corrections(self):
        return self._corrections

    @property
    def charge(self):
        return self.defect.charge
    
    @property
    def multiplicity(self):
        return self.defect.multiplicity
    
    @multiplicity.setter
    def multiplicity(self,multiplicity):
        self.defect._multiplicity = multiplicity
        return
    
    @property
    def name(self):
        name = self.defect.name 
        if self.label:
            name += f"({self.label})"
        return name

    @property
    def label(self):
        return self._label

    @property
    def delta_atoms(self):
        """
        Dictionary with Element as keys and particle difference between defect structure
        and bulk structure as values.
        """
        comp_defect = self.defect.defect_composition
        comp_bulk = self.bulk_structure.composition
        return get_delta_atoms_from_comp(comp_defect,comp_bulk)

    @property
    def classname(self):
        return self.__class__.__name__


    @staticmethod
    def from_jobs(job_defect, job_bulk, corrections, defect_structure=None,multiplicity=None,label=None):
        """
        Generate SingleDefectEntry object from VaspJob objects.

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
            Multiplicity of defect within the supercell. The default is None.
            If not provided is calculated by Pymatgen analysing the symmetry of the structure.
        label : (str), optional
            Additional label to add to defect specie. Does not influence non equilibrium calculations.
        Returns
        -------
        SingleDefectEntry
        """
        defect_structure = defect_structure if defect_structure else job_defect.initial_structure
        bulk_structure = job_bulk.final_structure
        energy_diff = job_defect.final_energy - job_bulk.final_energy
        charge = job_defect.charge
        
        return SingleDefectEntry.from_structures(defect_structure, bulk_structure, energy_diff, 
                                                 corrections,charge,multiplicity,label)


    @staticmethod
    def from_structures(defect_structure,bulk_structure,energy_diff,corrections,charge=0,multiplicity=None,label=None):
        """
        Generate SingleDefectEntry object from Structure objects.

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
            multiplicity of defect within the supercell. The default is None.
            If not provided is calculated by Pymatgen analysing the symmetry of the structure.
        label : (str), optional
            Additional label to add to defect specie. Does not influence non equilibrium calculations.

        Returns
        -------
        SingleDefectEntry
        """
        dsite,dtype = defect_finder(defect_structure, bulk_structure)
        module = importlib.import_module("pymatgen.analysis.defects.core")
        defect_class = getattr(module,dtype)
        defect = defect_class(bulk_structure,dsite,charge,multiplicity=multiplicity)
        
        return SingleDefectEntry(defect, bulk_structure, energy_diff, corrections,label)


    def __repr__(self):
        output = [
            "SingleDefectEntry :",
            "Defect = %s %s %s" %(self.defect.__class__.__name__,self.defect.site.specie.symbol,self.defect.site.frac_coords.__str__()),
            "Bulk System = %s" %self.bulk_structure.composition,
            "Energy = %.4f" %self.energy_diff,
            "Corrections = %.4f" %sum([v for v in self.corrections.values()]),
            "Charge = %i" %self.charge,
            "Multiplicity = %i" %self.multiplicity,
            "Name = %s" %self.name,
            "\n"
            ]
        return "\n".join(output)
    
    def __str__(self):
        return self.__repr__()


    def as_dict(self):
        """
        Returns:
            Json-serializable dict representation of SignleDefectData
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "defect": self.defect.as_dict(),
             "bulk_structure": self.bulk_structure.as_dict(),
             "energy_diff": self.energy_diff,
             "corrections": self.corrections,
             "label": self.label
             }
        return d        

    @classmethod
    def from_dict(cls,d):
        """
        Reconstitute a SingleDefectData object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of SingleDefectData.

        Returns:
            SingleDefectData object
        """
        defect = MontyDecoder().process_decoded(d['defect'])
        bulk_structure = Structure.from_dict(d['bulk_structure'])
        energy_diff = d['energy_diff']
        corrections = d['corrections']
        label = d['label'] if 'label' in d.keys() else None # ensure compatibility with old dictionaries
        return cls(defect,bulk_structure,energy_diff,corrections,label)


    def copy(self):
        return SingleDefectEntry(self.defect, self.bulk_structure, self.energy_diff, self.corrections,self.label)


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
    
    def defect_concentration(self, vbm, chemical_potentials, temperature=300, fermi_level=0.0):
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
        n = self.multiplicity * 1e24 / self.bulk_structure.volume
        conc = n * np.exp(-1.0 * self.formation_energy(vbm, chemical_potentials, fermi_level=fermi_level) /
                          (kb * temperature))
        return conc
    
    
    
class DefectComplexEntry:
    
    def __init__(self,defect_list,bulk_structure,energy_diff,corrections,charge=0,multiplicity=None,label=None):
        """
        Initializes the data of a defect complex. Inspired from the DefectEntry class in pymatgen (pmg).

        Args:
            defect_list: list Pymatgen defect objects (Vacancy, Interstitial or Substitution)
            bulk_structure: Pymatgen Structure without any defects
            energy_diff (float): difference btw energy of defect structure and energy of pure structure
            corrections (dict): Dict of corrections for defect formation energy. All values will be summed and
                                added to the defect formation energy.   
            charge : (int), optional
                Charge of the defect system. The default is 0.
            multiplicity : (int), optional
                multiplicity of defect within the supercell. The default is None.
                If not provided is calculated by Pymatgen analysing the symmetry of the structure.
            label : (str), optional
                Additional label to add to defect specie. Does not influence non equilibrium calculations.
        """        
        self._defect_list = defect_list 
        self._bulk_structure = bulk_structure
        self._energy_diff = energy_diff
        self._corrections = corrections if corrections else {}
        self._charge = int(charge)
        self._multiplicity = multiplicity if multiplicity else 1
        self._label = label 


    @property
    def defect_list(self):
        return self._defect_list
        
    @property
    def bulk_structure(self):
        return self._bulk_structure
    
    @property
    def energy_diff(self):
        return self._energy_diff
    
    @property
    def corrections(self):
        return self._corrections

    @property
    def charge(self):
        return self._charge
    
    @property
    def multiplicity(self):
        return self._multiplicity
    
    @multiplicity.setter
    def multiplicity(self,multiplicity):
        self._multiplicity = multiplicity
        return
    
    @property
    def name(self):
        name = '-'.join(self.defect_list_names)
        name = name + '_mult%i' %self.multiplicity
        if self.label:
            name += f"({self.label})"
        return name
    
    @property
    def label(self):
        return self._label
    
    @property
    def defect_list_names(self):
        names = []
        for d in self.defect_list:
            names.append(d.name.split('_mult', 1)[0])
        return names
    
    @property
    def delta_atoms(self):
        comp_bulk = self.bulk_structure.composition
        da_global = None
        for d in self.defect_list:
            da_single = get_delta_atoms_from_comp(d.defect_composition, comp_bulk)
            if da_global is None:
                da_global = da_single.copy()
            else:
                for e in da_single:
                    prec = da_global[e] if e in da_global.keys() else 0
                    da_global[e] = prec + da_single[e]
        
        return da_global


    @property
    def classname(self):
        return self.__class__.__name__
  
    
    @staticmethod
    def from_jobs(job_defect, job_bulk, corrections, defect_structure=None,multiplicity=None,label=None):
        """
        Generate DefectComplexEntry object from VaspJob objects.

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
            Multiplicity of defect within the supercell. The default is None.
            If not provided is calculated by Pymatgen analysing the symmetry of the structure.

        Returns
        -------
        DefectComplexEntry
        """ 
        defect_structure = defect_structure if defect_structure else job_defect.initial_structure
        bulk_structure = job_bulk.final_structure
        energy_diff = job_defect.final_energy - job_bulk.final_energy
        charge = job_defect.charge
        
        return DefectComplexEntry.from_structures(defect_structure, bulk_structure, energy_diff, corrections,
                                                  charge,multiplicity,label)
        
    
    @staticmethod
    def from_structures(defect_structure,bulk_structure,energy_diff,corrections,charge=0,multiplicity=None,label=None):
        """
        Generate DefectComplexEntry object from Structure objects.

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
            multiplicity of defect within the supercell. The default is None.
            If not provided is calculated by Pymatgen analysing the symmetry of the structure.

        Returns
        -------
        DefectComplexEntry
        """
        defect_list = []
        defects = defect_finder(defect_structure, bulk_structure)
        for dsite,dtype in defects:
            module = importlib.import_module("pynter.defects.pmg_defects_core") #modified pymatgen version with corrected multiplicity bug in super().init for Interstitials
            defect_class = getattr(module,dtype)
            defect = defect_class(bulk_structure,dsite,charge=charge,multiplicity=1)
            defect_list.append(defect)
        
        return DefectComplexEntry(defect_list, bulk_structure, energy_diff, corrections, charge, multiplicity,label)
    
    
    def __repr__(self):
        defect_list_string = [(d.__class__.__name__,d.site.specie.symbol,d.site.frac_coords.__str__()) for d in self.defect_list]
        output = [
            "DefectComplexEntry :",
            f"Defects =  {defect_list_string}",
            "Bulk System = %s" %self.bulk_structure.composition,
            "Energy = %.4f" %self.energy_diff,
            "Corrections = %.4f" %sum([v for v in self.corrections.values()]),
            "Charge = %i" %self.charge,
            f"Multiplicity = {self.multiplicity}",
            "Name = %s" %self.name,
            "\n"
            ]
        return "\n".join(output)
    
    def __str__(self):
        return self.__repr__()
    
    
    def as_dict(self):
        """
        Returns:
            Json-serializable dict representation of DefectComplexData
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "defect_list": [d.as_dict() for d in self.defect_list],
             "bulk_structure": self.bulk_structure.as_dict(),
             "energy_diff": self.energy_diff,
             "corrections": self.corrections,
             "charge": self.charge,
             "multiplicity": self.multiplicity,
             "label": self.label
             }
        return d        

    @classmethod
    def from_dict(cls,d):
        """
        Reconstitute a DefectComplexData object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of DefectComplexData.

        Returns:
            DefectComplexData object
        """
        defect_list = [MontyDecoder().process_decoded(d) for d in d['defect_list']]
        bulk_structure = Structure.from_dict(d['bulk_structure'])
        energy_diff = d['energy_diff']
        corrections = d['corrections']
        charge = d['charge']
        multiplicity = d['multiplicity']
        label = d['label'] if 'label' in d.keys() else None # ensure compatibility with old dictionaries
        return cls(defect_list,bulk_structure,energy_diff,corrections,charge,multiplicity,label)
    
    
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
    
    def defect_concentration(self, vbm, chemical_potentials, temperature=300, fermi_level=0.0):
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
        n = self.multiplicity * 1e24 / self.bulk_structure.volume
        conc = n * np.exp(-1.0 * self.formation_energy(vbm, chemical_potentials, fermi_level=fermi_level) /
                          (kb * temperature))
        return conc    
    
    
    