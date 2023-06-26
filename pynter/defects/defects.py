#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:08:42 2023

@author: villa
"""
from abc import ABCMeta, abstractmethod, abstractproperty
from monty.json import MSONable
from pymatgen.core.composition import Composition
from pymatgen.core.sites import PeriodicSite
import importlib
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pynter.tools.structure import is_site_in_structure, is_site_in_structure_coords, remove_oxidation_state_from_site
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
     
    @abstractproperty
    def defect_composition(self):
        """
        Defect composition as a Composition object
        """
        return

    @abstractproperty
    def defect_site_index(self):
        """
        Index of the defect site in the structure
        """
        return
    
    @property
    def defect_specie(self):
        return self.name.dspecie
    
    @property
    def defect_structure(self):
        """
        Structure of the defect
        """
        return self.generate_defect_structure() 
    
    @property
    def defect_type(self):
        return self.__class__.__name__
    
    @property
    def delta_atoms(self):
        """
        Dictionary with defect element symbol as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return

    @abstractmethod
    def get_multiplicity(self):
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
      
    @abstractproperty  
    def name(self):
        """
        Name of the defect
        """
        return

    @property
    def site(self):
        """
        Defect position as a site object
        """
        return self._defect_site
   
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
        Composition of the defect.
        """
        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] -= 1
        return Composition(temp_comp)
    
    @property
    def defect_site_index(self):
        """
        Index of the defect site in the bulk structure
        """
        _,index = is_site_in_structure_coords(self.site, self.bulk_structure)
        return index
    
    @property
    def delta_atoms(self):
        """
        Dictionary with element symbol as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return {self.site.specie.symbol:-1}
    
    def generate_defect_structure(self,bulk_structure=None,defect_site_index=None):
        """
        Generate a structure containing the defect starting from a bulk structure.

        Parameters
        ----------
        bulk_structure : (Structure), optional
            Bulk Structure. If not provided self.bulk_structure is used. The default is None.
        defect_site_index : (Structure), optional
            Index of the defect site in the bulk structure. If not provided 
            self.defect_site_index is used. The default is None.

        Returns
        -------
        structure : (Structure)
            Structure containing the defect.
        """
        bulk_structure = bulk_structure if bulk_structure else self.bulk_structure
        defect_site_index = defect_site_index if defect_site_index else self.defect_site_index
        defect_structure = bulk_structure.copy()
        defect_structure.remove_sites([defect_site_index])
        return defect_structure
    
    def get_multiplicity(self,**kwargs):
        """
        Get multiplicity of the defect in the structure

        Parameters
        ----------
        **kwargs : (dict)
            Arguments to pass to SpacegroupAnalyzer (symprec, angle_tolerance)
        """
        sga = SpacegroupAnalyzer(self.bulk_structure,**kwargs)
        symmetrized_structure = sga.get_symmetrized_structure()
        equivalent_sites = symmetrized_structure.find_equivalent_sites(self.site)
        return len(equivalent_sites)

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


    @property  
    def defect_composition(self):
        """
        Composition of the defect.
        """
        defect_index = self.defect_site_index
        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] += 1
        temp_comp[str(self.bulk_structure[defect_index].specie)] -= 1
        return Composition(temp_comp)

    @property
    def defect_site_index(self):
        """
        Index of the defect site in the bulk structure
        """
        return self.bulk_structure.index(self.site_in_bulk)

    @property
    def delta_atoms(self):
        """
        Dictionary with element symbol as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return {self.site.specie.symbol:1, self.site_in_bulk.specie.symbol:-1}
    
    def generate_defect_structure(self,bulk_structure=None,defect_site_index=None):
        """
        Generate a structure containing the defect starting from a bulk structure.

        Parameters
        ----------
        bulk_structure : (Structure), optional
            Bulk Structure. If not provided self.bulk_structure is used. The default is None.
        defect_site_index : (Structure), optional
            Index of the defect site in the bulk structure. If not provided 
            self.defect_site_index is used. The default is None.

        Returns
        -------
        structure : (Structure)
            Structure containing the defect.
        """
        bulk_structure = bulk_structure if bulk_structure else self.bulk_structure
        defect_site_index = defect_site_index if defect_site_index else self.defect_site_index
        defect_structure = bulk_structure.copy()
        defect_structure.replace(defect_site_index,self.defect_specie)  
        return defect_structure

    def get_multiplicity(self,**kwargs):
        """
        Get multiplicity of the defect in the structure

        Parameters
        ----------
        **kwargs : (dict)
            Arguments to pass to SpacegroupAnalyzer (symprec, angle_tolerance)
        """
        sga = SpacegroupAnalyzer(self.bulk_structure,**kwargs)
        symmetrized_structure = sga.get_symmetrized_structure()
        equivalent_sites = symmetrized_structure.find_equivalent_sites(self.site_in_bulk)
        return len(equivalent_sites)

    @property 
    def name(self):
        """
        Name for this defect. Behaves like a string with additional attributes.
        """
        return DefectName(self.defect_type,self.site.specie.symbol,self.site_in_bulk.specie.symbol,self.label)   
    
    @property
    def site_in_bulk(self):
        return self._site_in_bulk
    
    def get_site_in_bulk(self):
        try:
            site = min(
                self.bulk_structure.get_sites_in_sphere(self.site.coords, 0.5, include_index=True),
                           key=lambda x: x[1])  
            # there's a bug in pymatgen PeriodicNeighbour.from_dict and the specie attribute, get PeriodicSite instead
            site = PeriodicSite(site.species, site.frac_coords, site.lattice)
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
        Composition of the defect.
        """
        temp_comp = self.bulk_structure.composition.as_dict()
        temp_comp[str(self.site.specie)] += 1
        return Composition(temp_comp)
    
    @property
    def defect_site_index(self):
        """
        Index of the defect site in the defect structure
        """
        is_site,index = is_site_in_structure(self.site, self.defect_structure)
        return index  # more flexibility with is_site_in_structure
       
    @property
    def delta_atoms(self):
        """
        Dictionary with element symbol as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return {self.site.specie.symbol:1}
    
    def generate_defect_structure(self,bulk_structure=None):
        """
        Generate a structure containing the defect starting from a bulk structure.

        Parameters
        ----------
        bulk_structure : (Structure), optional
            Bulk Structure. If not provided self.bulk_structure is used. The default is None.

        Returns
        -------
        structure : (Structure)
            Structure containing the defect.
        """
        bulk_structure = bulk_structure if bulk_structure else self.bulk_structure
        defect_structure = bulk_structure.copy()
        defect_structure.append(self.site.species,self.site.frac_coords)
        return defect_structure
        

    def get_multiplicity(self):
        raise NotImplementedError('Not implemented for Interstitial')

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
    
    def __init__(self,defect_site,bulk_structure,charge=None,multiplicity=None,
                 label=None,defect_structure=None):
        """
        defect_structure: (Structure)
            Structure containing the polaron. If not provided the site index is searched 
            in the bulk structure, and the defect_structure is set equal to the bulk structure.
        """
        super().__init__(defect_site,bulk_structure,charge,multiplicity,label) 
        self._defect_structure = defect_structure
        
        
    @property
    def defect_composition(self):
        """
        Composition of the defect
        """
        return self.bulk_structure.composition
    
    @property
    def defect_site_index(self):
        """
        Index of the defect site in the structure
        """
        return self.defect_structure.index(self.site)
        
    @property
    def delta_atoms(self):
        """
        Dictionary with delement as keys and difference in particle number 
        between defect and bulk structure as values
        """
        return {}
    
    def generate_defect_structure(self,bulk_structure=None):
        """
        Structure containing the polaron. If not provided the site index is searched 
        in the bulk structure, and the defect_structure is set equal to the bulk structure.
        """
        if self._defect_structure:
            return self._defect_structure
        else:
            bulk_structure = bulk_structure if bulk_structure else self.bulk_structure
            return bulk_structure  
        
        
    def get_multiplicity(self,**kwargs):
        """
        Get multiplicity of the defect in the structure

        Parameters
        ----------
        **kwargs : (dict)
            Arguments to pass to SpacegroupAnalyzer (symprec, angle_tolerance)
        """
        sga = SpacegroupAnalyzer(self.bulk_structure,**kwargs)
        symmetrized_structure = sga.get_symmetrized_structure()
        equivalent_sites = symmetrized_structure.find_equivalent_sites(self.site)
        return len(equivalent_sites)

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
    def defect_composition(self):
        return self.defect_structure.composition
        
    @property
    def defect_names(self):
        return [d.name for d in self.defects]

    @property
    def defect_structure(self):
        return self.generate_defect_structure()

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

    def generate_defect_structure(self,bulk_structure=None):
        """
        Generate a structure containing the defect starting from a bulk structure.
        If not provided self.bulk_structure is used.
        """
        bulk_structure = bulk_structure if bulk_structure else self.bulk_structure
        structure = bulk_structure.copy()
        for df in self.defects:
            df_structure = df.generate_defect_structure(structure)
            structure = df_structure.copy()
        return structure

    def get_multiplicity(self):
        raise NotImplementedError('Not implemented for DefectComplex')

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
        return self.fullname.__str__()
    
    def __repr__(self):
        return self.fullname.__repr__()
    
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
        return self.fullname.__str__()
    
    def __repr__(self):
        return self.fullname.__repr__()
    
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


def get_defect_name_from_string(string):
    if '-' in string:
        return DefectComplexName.from_string(string)
    else:
        return DefectName.from_string(string)


def create_interstitials(structure,elements,supercell_size=None,**kwargs):
    """
    Create Interstitial objects based on Voronoi with pymatgen. 

    Parameters
    ----------
    structure : (Structure)
        Bulk structure.
    elements : (list)
        List of element symbols.
    supercell_size : (int), optional
        Input for the make_supercell function of the Structure class.
        If None the input structure is not modified. The default is None.
    kwargs: 
        Arguments to pass to VoronoiInterstitialGenerator:
            clustering_tol: Tolerance for clustering the Voronoi nodes.
            min_dist: Minimum distance between an interstitial and the nearest atom.
            ltol: Tolerance for lattice matching.
            stol: Tolerance for structure matching.
            angle_tol: Angle tolerance for structure matching.
            kwargs: Additional keyword arguments for the ``TopographyAnalyzer`` constructor.

    Returns
    -------
    defects : (list)
        List of Interstitial objects
    """
    from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
    
    defects = []
    bulk_structure = structure.copy()
    if supercell_size:
        bulk_structure.make_supercell(supercell_size)
    generator = VoronoiInterstitialGenerator().generate(bulk_structure,elements)
    for inter in generator:
        bulk_structure.remove_oxidation_states()
        remove_oxidation_state_from_site(inter.site)
        interstitial = Interstitial(inter.site, bulk_structure,multiplicity=inter.multiplicity,
                                    label=f'mult{inter.multiplicity}')
        defects.append(interstitial)
    return defects


def create_substitutions(structure,elements_to_replace,supercell_size=None,site_indexes=None):
    """
    Create Substitution objects starting from a bulk structure (unit cell or supercell).

    Parameters
    ----------
    structure : Structure
        Bulk structure, both unit cell or supercell can be used as input.
    elements_to_replace : (str), optional
        Dict with element symbol of specie to be replaced as keys and element 
        symbol of the species to be replaced with as values ({'old_El':'new_El'}).
    supercell_size : (int or numpy array), optional
        Input for the generate_defect_structure function of the old pymatgen Substitution class.
    site_indexes : (list), optional
        Site indexes where the defects need to be created. If None the first site found is 
        used. The default is None.

    Returns
    -------
    defects : (list)
        List of Substitution objects
    """
    defects = [] 
    bulk_structure = structure.copy()
    if supercell_size:
        bulk_structure.make_supercell(supercell_size)
    if site_indexes:
        sites = [bulk_structure[i] for i in site_indexes]
    else:
        sites = bulk_structure.sites
    for el_to_sub,el_subbed in elements_to_replace.items():
        for site in sites:
            if site.specie.symbol == el_to_sub:
                defect_site = PeriodicSite(el_subbed,site.frac_coords,site.lattice)
                defects.append(Substitution(defect_site, bulk_structure,site_in_bulk=site))
                if not site_indexes:
                    break
    return defects   


def create_vacancies(structure,elements=None,supercell_size=None,site_indexes=None):
    """
    Create structures with vacancies starting from a bulk structure (unit cell or supercell).

    Parameters
    ----------
    structure : Structure
        Bulk structure, both unit cell or supercell can be used as input.
    elements : (str), optional
        Symbol of the elements for which vacancies are needed.
        If None all of the elements are considered. The default is None.
    supercell_size : (int or numpy array), optional
        Input for the make_supercell function of the Structure class.
        If None the input structure is not modified. The default is None.
    site_indexes : (list), optional
        Site indexes where the defects need to be created. If None the first site found is 
        used. The default is None.

    Returns
    -------
    defects : (list)
        List of Vacancy objects
    """
    defects = []
    bulk_structure = structure.copy()
    if supercell_size:
        bulk_structure.make_supercell(supercell_size)
    if not elements:
        elements = [el.symbol for el in bulk_structure.composition.elements]
    if site_indexes:
        sites = [bulk_structure[i] for i in site_indexes]
    else:
        sites = bulk_structure.sites
    for el in bulk_structure.composition.elements:
        for site in sites:
            if el.symbol in elements:
                if site.specie == el:
                    vacancy = Vacancy(site, bulk_structure)
                    defects.append(vacancy)
                    if not site_indexes:
                        break
    return defects        


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
    mod_label = label + '$'
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
    

def get_old_pmg_object(defect):
    """
    Convert defect object to the old version of pymatgen defect object.
    """
    module = importlib.import_module("pynter.defects.pmg.pmg_defects_core")
    pmg_class = getattr(module,defect.defect_type)
    charge = defect.charge if defect.charge else 0
    pmg_object = pmg_class(structure=defect.bulk_structure,defect_site=defect.site,
                           charge=charge,multiplicity=defect.multiplicity)
    return pmg_object

def get_pmg_object(defect):
    """
    Convert defect object to pymatgen defect object.
    """
    module = importlib.import_module("pymatgen.analysis.defects.core")
    pmg_class = getattr(module,defect.defect_type)
    pmg_object = pmg_class(structure=defect.bulk_structure,site=defect.site,
                           multiplicity=defect.multiplicity)
    return pmg_object
    