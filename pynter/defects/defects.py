#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:08:42 2023

@author: villa
"""
from abc import ABCMeta, abstractmethod, abstractproperty
from monty.json import MSONable
from pymatgen.core.composition import Composition
import json
from monty.json import MontyDecoder, MontyEncoder
from pymatgen.core.sites import PeriodicSite

# Adapted from pymatgen.analysis.defect.core 


class Defect(MSONable,metaclass=ABCMeta): #MSONable contains as_dict and from_dict methods
    """
    Abstract class for a single point defect
    """

    def __init__(self, defect_site, bulk_structure, charge=None, multiplicity=None,label=None):
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
        self._charge = charge
        self._multiplicity = multiplicity 
        self._label = label


    def __repr__(self):
        return "{} : {} {}".format(self.defect_type,self.site.frac_coords, self.name.dspecie)
    
    def __print__(self):
        return self.__repr__()

    @property
    def bulk_structure(self):
        """
        Structure without defects.
        """
        return self._bulk_structure
    
    @property
    def charge(self):
        """
        Charge of the defect.
        """
        return self._charge
    
    @property 
    @abstractmethod
    def defect_composition(self):
        """
        Defect composition as a Composition object
        """
        return

    @property
    def site(self):
        """
        Defect position as a site object
        """
        return self._defect_site

    @property
    def defect_specie(self):
        return self.name.dspecie
    
    @property
    def defect_type(self):
        return self.__class__.__name__
    
    @property
    @abstractmethod
    def delta_atoms(self):
        """
        Dictionary with defect element as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return

    @property
    def label(self):
        return self._label

    @property
    def multiplicity(self):
        """
        Multiplicity of a defect site within the structure
        """
        return self._multiplicity
      
    @property  
    @abstractmethod
    def name(self):
        """
        Name of the defect
        """
        return
   
    @property
    def site_concentration_in_cm3(self):
        return self.multiplicity * 1e24 / self.bulk_structure.volume 
   
    @property  
    def symbol(self):
        """
        Latex formatted name of the defect
        """
        return self.name.symbol

    @property
    def symbol_with_charge(self):
        """
        Name in latex format with charge written as a number.
        """
        return format_legend_with_charge_number(self.symbol,self.charge)
    
    @property
    def symbol_with_charge_kv(self):
        """
        Name in latex format with charge written with kroger and vink notation.
        """
        return format_legend_with_charge_kv(self.symbol,self.charge)

    def set_charge(self, new_charge=0.0):
        """
        Sets the charge of the defect.
        """
        self._charge = new_charge
        return
    
    def set_label(self,new_label):
        """
        Sets the label of the defect
        """
        self._label = new_label
        return

    def set_multiplicity(self, new_multiplicity=1):
        """
        Sets the charge of the defect.
        """
        self._multiplicity = new_multiplicity
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
    def delta_atoms(self):
        """
        Dictionary with delement as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return {self.site.specie:-1}

    @property
    def name(self):
        """
        Name of the defect. Behaves like a string with additional attributes.
        """
        return DefectName(self.defect_type,self.site.specie.symbol,None,self.label)
        
    
class Substitution(Defect):
    """
    Subclass of Defect for substitutional defects.
    """

    def __init__(self,defect_site,bulk_structure,charge=None,multiplicity=None,label=None,site_in_bulk=None):
        """
        site_in_bulk: (PeriodicSite)
            Original Site in bulk structure were substitution took place
        """
        super().__init__(defect_site,bulk_structure,charge,multiplicity,label)  
        self._site_in_bulk = site_in_bulk if site_in_bulk else self.get_site_in_bulk()


    @property  # type: ignore
    def defect_composition(self):
        """
        Returns: Composition of defect.
        """
        defect_index = self.site_in_bulk.index
        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] += 1
        temp_comp[str(self.bulk_structure[defect_index].specie)] -= 1
        return Composition(temp_comp)

    @property
    def delta_atoms(self):
        """
        Dictionary with delement as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return {self.site.specie:1, self.site_in_bulk.specie:-1}

    @property 
    def name(self):
        """
        Name for this defect. Behaves like a string with additional attributes.
        """
        return DefectName(self.defect_type,self.site.specie.symbol,self.site_in_bulk.specie,self.label)   
    
    @property
    def site_in_bulk(self):
        return self._site_in_bulk
    
    def get_site_in_bulk(self):
        try:
            site = min(
                self.bulk_structure.get_sites_in_sphere(self.site.coords, 0.5, include_index=True),
                           key=lambda x: x[1])  
            # there's a bug in pymatgen PeriodicNeighbour._from_dict and the specie attribute, get PeriodicSite instead
            site = PeriodicSite(site.species, site.coords, site.lattice)
            return site
        except:
            return ValueError("""No equivalent site has been found in bulk, defect and bulk structures are too different.\
Try using the unrelaxed defect structure or provide bulk site manually""")
    
 
    
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
    def delta_atoms(self):
        """
        Dictionary with delement as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return {self.site.specie:1}

    @property
    def name(self):
        """
        Name of the defect. Behaves like a string with additional attributes.
        """
        return DefectName(self.defect_type,self.site.specie.symbol,None,self.label)
        
      
class Polaron(Defect):
    """
    Subclass of Defect for polarons
    """    
    @property
    def defect_composition(self):
        return self.bulk_structure.composition
    
    @property
    def delta_atoms(self):
        """
        Dictionary with delement as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return {}

    @property
    def name(self):
        """
        Name of the defect. Behaves like a string with additional attributes.
        """
        return DefectName(self.defect_type,self.site.specie.symbol,None,self.label)
  
    
    
class DefectComplex(MSONable,metaclass=ABCMeta):

    def __init__(self,defects,bulk_structure,charge=None,multiplicity=None,label=None):
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
        self._label = label


    def __repr__(self):
        return "{} : {}".format(self.defect_type,
        [(d.defect_type,d.name.dspecie,d.site.frac_coords.__str__()) 
         for d in self.defects])

    def __str__(self):
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
    def defects(self):
        """
        List of single defects consituting the complex.
        """
        return self._defects
        
    @property
    def defect_names(self):
        return [d.name for d in self.defects]

    @property
    def defect_type(self):
        return self.__class__.__name__
    
    @property
    def delta_atoms(self):
        """
        Dictionary with Element as keys and particle difference between defect structure
        and bulk structure as values.
        """
        da_global = None
        for d in self.defects:
            da_single = d.delta_atoms
            if da_global is None:
                da_global = da_single.copy()
            else:
                for e in da_single:
                    prec = da_global[e] if e in da_global.keys() else 0
                    da_global[e] = prec + da_single[e]       
        return da_global 

    @property
    def label(self):
        return self._label

    @property
    def multiplicity(self):
        """
        Multiplicity of a defect site within the structure
        """
        return self._multiplicity

    @property
    def name(self):
        """
        Name of the defect. Behaves like a string with additional attributes.
        """
        return DefectComplexName([d.name for d in self.defects],self.label)

    @property
    def site_concentration_in_cm3(self):
        return self.multiplicity * 1e24 / self.bulk_structure.volume 
    
    @property
    def sites(self):
        return [d.site for d in self.defects]
    
    @property
    def symbol(self):
        return self.name.symbol

    @property
    def symbol_with_charge(self):
        """
        Name in latex format with charge written as a number.
        """
        return format_legend_with_charge_number(self.symbol,self.charge)
    
    @property
    def symbol_with_charge_kv(self):
        """
        Name in latex format with charge written with kroger and vink notation.
        """
        return format_legend_with_charge_kv(self.symbol,self.charge)


    def set_charge(self, new_charge):
        """
        Sets the charge of the defect.
        """
        self._charge = new_charge
        return    

    def set_label(self,new_label):
        """
        Sets the label of the defect
        """
        self._label = new_label
        return

    def set_multiplicity(self, new_multiplicity):
        """
        Sets the charge of the defect.
        """
        self._multiplicity = new_multiplicity
        return
    

class DefectName(str,MSONable):
    
    def __new__(cls, dtype,dspecie,bulk_specie=None,label=None):
        if dtype == 'Vacancy':
            name = 'Vac_%s' %dspecie
            symbol = '$V_{%s}$' %dspecie
            
        elif dtype == 'Substitution':
            name = 'Sub_%s_on_%s' %(dspecie,bulk_specie)
            symbol = '$%s_{%s}$' %(dspecie,bulk_specie)
            
        elif dtype == 'Interstitial':
            name = 'Int_%s' %dspecie
            symbol = '$%s_{i}$' %dspecie
            
        elif dtype == 'Polaron':
            name = 'Pol_%s' %dspecie
            symbol = '$%s_{%s}$' %(dspecie,dspecie)
            
        if label:
            fullname = name +'(%s)'%label
            symbol = symbol +'(%s)'%label
        else:
            fullname = name
        
        #get str instance
        instance = super().__new__(cls,fullname) 
        instance._dtype = dtype
        instance._dspecie = dspecie
        instance._bulk_specie = bulk_specie
        instance._name = name
        instance._label = label
        instance._fullname = fullname
        instance._symbol = symbol
        
        return instance


    def __init__(self,dtype,dspecie,bulk_specie=None,label=None):
        """
        Class to systematically organize the name of a defect.
        Behaves like a string, with the additional attributes.

        Parameters
        ----------
        dtype : (str)
            Type of the defect. Types are "Vacancy","Substitution","Interstitial","Polaron".
        dspecie : (str)
            Element of defect specie.
        bulk_specie : (str)
            Element of the specie that has been replaced. To be set only for "Substitution".
        label : (str), optional
            Label for the defect. Will be displayed in parenthesis after the name. The default is None.
        """
        super().__init__() # init string


    def __str__(self):
        return self.fullname
    
    def __repr__(self):
        return self.fullname
    
    def __iter__(self):
        return [self].__iter__() # dummy iter to handle single defecs and complexes the same way 

    def __eq__(self, other):
        if isinstance(other, str):
            return self.fullname == other
        elif isinstance(other, DefectName):
            return self.fullname == other.fullname
        else:
            return False

    def __hash__(self):
        return hash(self.fullname)

    @property    
    def bulk_specie(self):
        if hasattr(self, '_bulk_specie'):
            bulk_specie = self._bulk_specie
        else:
            bulk_specie = None
        return bulk_specie

    @property        
    def dspecie(self):
        return self._dspecie

    @property
    def dtype(self):
        return self._dtype

    @property
    def fullname(self):
        return self._fullname

    @property
    def label(self):
        return self._label

    @property
    def name(self):
        return self._name

    @property
    def symbol(self):
        return self._symbol
    
    
    @staticmethod
    def from_string(string):
        args = DefectName._get_args_from_string(string)
        return DefectName(*args)
          
    def _get_args_from_string(string):
        if '(' in string:
            name,label = string.split('(')
            label = label.strip(')')
        else:
            name = string
            label=None
        nsplit = name.split('_')
        ntype = nsplit[0]
        el = nsplit[1]
        bulk_specie = None
        if ntype=='Vac':
            dtype = 'Vacancy'
            dspecie = el
        elif ntype=='Int':
            dtype = 'Interstitial'
            dspecie = el
        elif ntype=='Sub':
            dtype = 'Substitution'
            dspecie = el
            el_bulk = nsplit[3]
            bulk_specie = el_bulk
        elif ntype=='Pol':
            dtype = 'Polaron'
            dspecie = el
        return dtype,dspecie,bulk_specie,label
            
            
class DefectComplexName(str,MSONable):
    
    def __new__(cls, defect_names, label=None):
        name = '-'.join([n.name for n in defect_names]) #exclude labels
        if label:
            fullname =  name + '(%s)'%label
        else:
            fullname = name
        instance = super().__new__(cls,fullname)        
        instance._name = name
        instance._label = label
        instance._fullname = fullname
        instance._defect_names = defect_names
        return instance
    
    def __init__(self,defect_names,label=None):
        """
        Class to systematically organize the name of a defect complex. The name of the 
        complex concatenates with '-' the names of the single defects.
        Behaves like a string, with the additional attributes.

        Parameters
        ----------
        defect_names : (list)
            List of DefectName objects.
        label : (str), optional
            Label for the defect. Will be displayed in parenthesis after the name. The default is None.
        """
        super().__init__()
        
                
    def __str__(self):
        return self.fullname
    
    def __repr__(self):
        return self.fullname
    
    def __iter__(self):
        return self.defect_names.__iter__()

    def __eq__(self, other):
        if isinstance(other, str):
            return self.fullname == other
        elif isinstance(other, DefectName):
            return self.fullname == other.fullname
        else:
            return False

    def __hash__(self):
        return hash(self.fullname)

    @property
    def defect_names(self):
        return self._defect_names

    @property
    def dspecies(self):
        return [n.dspecie for n in self.defect_names]
    
    @property
    def dtype(self):
        return 'DefectComplex'
        
    @property
    def fullname(self):
        return self._fullname

    @property
    def label(self):
        return self._label

    @property
    def name(self):
        return self._name
    
    @property
    def symbol(self):
        symbol = '-'.join([n.symbol.split('(')[0] for n in self.defect_names]) 
        if self.label:
            symbol = symbol + '(%s)'%self.label 
        return symbol
    
    
    @staticmethod
    def from_string(string):
        args = DefectComplexName._get_args_from_string(string)
        return DefectComplexName(*args)
    
    def _get_args_from_string(string):
        if '(' in string:
            name,label = string.split('(')
            label = label.strip(')')
        else:
            name = string
            label=None
        names = name.split('-')
        defect_names = [DefectName.from_string(n) for n in names]
        return defect_names,label
        

def format_legend_with_charge_number(label,charge):
    """
    Get label in latex format with charge written as a number.
   
    Parameters
    ----------
    label : (str)
        Original name of the defect.
    charge : (int or float)
        Charge of the defect. Floats are converted to integers.
    """
    s = label
    charge = int(charge)
    if charge > 0:
        q = '+'+str(charge)
    elif charge == 0:
        q = '\;' + str(charge)
    else:
        q = str(charge)
    return s + '$^{' + q + '}$'


def format_legend_with_charge_kv(label,charge):
    """
    Get label in latex format with charge written with kroger and vink notation.

    Parameters
    ----------
    label : (str)
        Original name of the defect.
    charge : (int or float)
        Charge of the defect. Floats are converted to integer.
    """
    mod_label = label[:-1]
    charge = int(charge)
    if charge < 0:
        for i in range(0,abs(charge)):
            if i == 0:
                mod_label = mod_label + "^{"
            mod_label = mod_label + "'"
        mod_label = mod_label + "}"
    elif charge == 0:
        mod_label = mod_label + "^{x}"
    elif charge > 0:
        for i in range(0,charge):
            if i == 0:
                mod_label = mod_label + "^{"
            mod_label = mod_label + "°"
        mod_label = mod_label + "}"
    
    mod_label = mod_label + "$"
    
    return mod_label

    
def get_delta_atoms(structure_defect,structure_bulk):
    """
    Function to build delta_atoms dictionary starting from Pymatgen Structure objects.
    ----------
    structure_defect : (Pymatgen Structure object)
        Defect structure.
    structure_bulk : (Pymatgen Structure object)
        Bulk structure.
    Returns
    -------
    delta_atoms : (dict)
        Dictionary with Element as keys and delta n as values.
    """
    comp_defect = structure_defect.composition
    comp_bulk = structure_bulk.composition
        
    return get_delta_atoms_from_comp(comp_defect, comp_bulk)


def get_delta_atoms_from_comp(comp_defect,comp_bulk):
    """
    Function to biuld delta_atoms dictionary starting from Pymatgen Structure objects.
    ----------
    comp_defect : (Pymatgen Composition object)
        Defect structure.
    comp_bulk : (Pymatgen Composition object)
        Bulk structure.
    Returns
    -------
    delta_atoms : (dict)
        Dictionary with Element as keys and delta n as values.
    """
    delta_atoms = {}
    for el,n in comp_defect.items():
        nsites_defect = n
        nsites_bulk = comp_bulk[el] if el in comp_bulk.keys() else 0
        delta_n = nsites_defect - nsites_bulk
        if delta_n != 0:
            delta_atoms[el] = delta_n
        
    return delta_atoms    
    
    