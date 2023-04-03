#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:08:42 2023

@author: villa
"""
from abc import ABCMeta, abstractmethod
from pymatgen.core.composition import Composition

# Adapted from pymatgen.analysis.defect.core 

class Defect:
    """
    Abstract class for a single point defect
    """

    def __init__(self, defect_site, bulk_structure, charge, multiplicity):
        """
        Base class for defect objets

        Parameters
        ----------
        defect_site : (Site)
            Pymatgen Site object of the defect.
        bulk_structure : (Structure)
            Pymatgen Structure without defects.
        charge : (int or float)
            Charge of the defect.
        multiplicity : (int)
            Multiplicity of defect within the simulation cell.
        """
        self._defect_site = defect_site
        self._bulk_structure = bulk_structure
        self._charge = int(charge)
        self._multiplicity = multiplicity 

    def __repr__(self):
        return "{} : {}".format(self.__class__.__name__,self.site)
    
    def __print__(self):
        return self.__repr__()

    @property
    def bulk_structure(self):
        """
        Returns the structure without any defects.
        """
        return self._bulk_structure

    @property
    def charge(self):
        """
        Charge of the defect.
        """
        return self._charge

    @property
    def site(self):
        """
        Defect position as a site object
        """
        return self._defect_site

    @property
    def multiplicity(self):
        """
        Multiplicity of a defect site within the structure
        """
        return self._multiplicity
      
    @property 
    @abstractmethod
    def defect_composition(self):
        """
        Defect composition as a Composition object
        """
        return

    @property  
    @abstractmethod
    def name(self):
        """
        Name of the defect
        """
        return

    def set_charge(self, new_charge=0.0):
        """
        Sets the charge of the defect.
        """
        self._charge = int(new_charge)
        return
        
        
class Vacancy(Defect):
    
    """
    Subclass of Defect for single vacancies.
    """

    @property
    def defect_composition(self):
        """
        Returns: Composition of defect.
        """
        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] -= 1
        return Composition(temp_comp)
    

    @property
    def name(self):
        """
        Name of the defect.
        """
        return "Vac_{}".format(self.site.specie)
    
    
    
class Substitution(Defect):
    """
    Subclass of Defect for substitutional defects.
    """

    @property  # type: ignore
    def defect_composition(self):
        """
        Returns: Composition of defect.
        """
        poss_deflist = min(
            self.bulk_structure.get_sites_in_sphere(self.site.coords, 0.1, include_index=True),
            key=lambda x: x[1],
        )
        defindex = poss_deflist[0][2]

        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] += 1
        temp_comp[str(self.bulk_structure[defindex].specie)] -= 1
        return Composition(temp_comp)


    @property 
    def name(self):
        """
        Returns a name for this defect
        """
        poss_deflist = min(
            self.bulk_structure.get_sites_in_sphere(self.site.coords, 0.1, include_index=True),
            key=lambda x: x[1],
        )
        defindex = poss_deflist[0][2]
        return "Sub_{}_on_{}".format(self.site.specie, self.bulk_structure[defindex].specie)    
    
    
    
class Interstitial(Defect):
    """
    Subclass of Defect for interstitial defects.
    """
    
    @property
    def defect_composition(self):
        """
        Returns: Defect composition.
        """
        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] += 1
        return Composition(temp_comp)
    
    @property
    def name(self):
        """
        Name of the defect
        """
        return "Int_{}".format(self.site.specie)    
    
    
class Polaron(Defect):
    """
    Subclass of Defect for polarons
    """    
    @property
    def defect_composition(self):
        return self.bulk_structure.composition
    
    @property
    def name(self):
        """
        Name of the defect
        """
        return "Pol_{}".format(self.site.specie)   
  
    
    
class DefectComplex:

    def __init__(self,defects,bulk_structure,charge,multiplicity):
        """
        Class to describe defect complexes

        Parameters
        ----------
        defects : (list)
            List of Defect objects.
        bulk_structure : (Structure)
            Pymatgen Structure of the bulk material.
        charge : (int or float)
            Charge of the defect.
        multiplicity : (int)
            Multiplicity of the defect.
        """
        self._defects = defects
        self._bulk_structure = bulk_structure
        self._charge = charge
        self._multiplicity = multiplicity

    @property
    def defects(self):
        """
        List of single defects consituting the complex.
        """
        return self._defects
        
    @property
    def bulk_structure(self):
        """
        Returns the structure without any defects.
        """
        return self._bulk_structure

    @property
    def charge(self):
        """
        Charge of the defect.
        """
        return self._charge

    @property
    def multiplicity(self):
        """
        Multiplicity of a defect site within the structure
        """
        return self._multiplicity
    
    
    
    
    
    