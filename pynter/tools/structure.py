#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 10:52:08 2020

@author: villa
"""

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view


def _is_site_in_structure_old(site,structure,tol=1e-03):
    """
    Check if Site is part of the Structure list. This function is needed because 
    sometimes doing a simple check ("site in structure") doesn't work. This function performes
    a check on the coordinates and the element on the site. Therefore it is more reliable.

    Parameters
    ----------
    site : (Site)
        PeriodicSite or Site object.
    structure : (Structure)
        Pymatgen Structure object.
    tol : (float), optional
        Tolerance for fractional coordinates. The default is 1e-03.

    Returns
    -------
    is_site_in_structure : (bool)
    index : (int)
        Index of site in structure in case site is_site_in_structure returns True
        If False index will be None.
    """
    is_site_in_structure = False
    for s in structure:
        if np.allclose(site.frac_coords,s.frac_coords,atol=tol) and site.specie.symbol == s.specie.symbol:
            is_site_in_structure = True
            index = structure.index(s)
            return is_site_in_structure,index
    index=None
    return is_site_in_structure,index


def is_site_in_structure(site,structure,tol=1e-03):
    """
    Check if Site is part of the Structure list. This function is needed because 
    sometimes doing a simple check ("site in structure") doesn't work. This function performes
    a check on the coordinates and the element on the site. Therefore it is more reliable.

    Parameters
    ----------
    site : (Site)
        PeriodicSite or Site object.
    structure : (Structure)
        Pymatgen Structure object.
    tol : (float), optional
        Tolerance for fractional coordinates. The default is 1e-03.

    Returns
    -------
    is_site_in_structure : (bool)
    index : (int)
        Index of site in structure in case site is_site_in_structure returns True
        If False index will be None.
    """
    is_site_in_structure,index = False,None
    check,index = is_site_in_structure_coords(site, structure,tol=tol)
    if check:
        s = structure[index]
        if site.specie.symbol == s.specie.symbol:
            is_site_in_structure =True
    
    return is_site_in_structure,index


def _is_site_in_structure_coords_old(site,structure,tol=1e-03):
    """
    Check if Site coordinates are prensent in the Structure list. 

    Parameters
    ----------
    site : (Site)
        PeriodicSite or Site object.
    structure : (Structure)
        Pymatgen Structure object.
    tol : (float), optional
        Tolerance for fractional coordinates. The default is 1e-03.

    Returns
    -------
    is_site_in_structure_coords : (bool)
    index : (int)
        Index of site in structure in case site is_site_in_structure_coords returns True
        If False index will be None.
    """
    is_site_in_structure_coords = False
    for s in structure:
        if np.allclose(site.frac_coords,s.frac_coords,atol=tol):
            is_site_in_structure_coords = True
            index = structure.index(s)
            return is_site_in_structure_coords,index
    index=None
    return is_site_in_structure_coords,index


def is_site_in_structure_coords(site,structure,tol=1e-03):
    """
    Check if Site coordinates are prensent in the Structure. Calculates distance between target site and 
    each site in reference structure, takes the minimum value and returns True if is within a tolerance.
    The periodicity is accounted for (considering the lattice associated to the PeriodicSite.
    The tolerance is normalized with respect to lattice vector size. 

    Parameters
    ----------
    site : (Site)
        PeriodicSite or Site object.
    structure : (Structure)
        Pymatgen Structure object.
    tol : (float), optional
        Tolerance for site comparison. The distance between sites in target and reference stucture is used, 
        periodicity is accounted for. The tolerance is normalized with respect to lattice vector size. 
        The default is 1e-03.

    Returns
    -------
    is_site_in_structure_coords : (bool)
    index : (int)
        Index of site in structure in case site is_site_in_structure_coords returns True
        If False index will be None.
    """
    l = site.lattice
    tol = np.sqrt(l.a**2 + l.b**2 + l.c**2) * tol #input is normalized with respect to lattice vector
    distances=[]
    for s in structure:
        distances.append(s.distance(site))
    
    distance_min = min(distances)
    if distance_min < tol:
        i_min = distances.index(distance_min)
        index = structure.index(structure[i_min])
        return True, index
    else:
        return False,None


def sort_sites_to_ref_coords(structure,structure_ref,extra_sites=[],tol=1e-03,get_indexes=False):
    """
    Sort Sites of one structure to match the order of coordinates in a reference structure. 

    Parameters
    ----------
    structure : (Structure)
        Structure to sort.
    structure_ref : (Structure)
        Reference Structure.
    extra_sites : (list), optional
        Sites to append at the end of the structure. The default is [].
    tol: (float)
        Tolerance for site comparison. The distance between sites in target and reference stucture is used, 
        periodicity is accounted for. The tolerance is normalized with respect to lattice vector size. 
        The default is 1e-03.
    get_indexes : (bool), optional
        Get list of mapping indexes for target structure sites in reference structure. The default is False.

    Returns
    -------
    new_structure : (Structure)
        Sorted Structure
    indexes : (list)
        If get_indexes is True a tuple is returned. List of mapping indexes for target structure sites
        in reference structure.
    """
    df = structure
    bk = structure_ref
    indexes = []
    new_sites=[]
    for s in df:
        check,index = is_site_in_structure_coords(s,bk,tol=tol)
        if check:
            indexes.append(index)
    for w in range(0,len(bk)):
        new_sites.insert(indexes[w],df[w])
    for s in extra_sites:
        new_sites.append(s)
            
    new_structure = df.copy()
    new_structure._sites = new_sites 
    if get_indexes:
        return new_structure, indexes
    else:
        return new_structure


def view_structure_with_ase(structures):
    """
    Visualize a Structure object or a list of Structure objects with the ASE. 
    First the Structure objects are converted into an ase Atom object, then "view" is used to visualize them.
    """
    if type(structures) == list:
        atoms=[]
        for s in structures:
            atoms.append(AseAtomsAdaptor.get_atoms(s))
    else:
        atoms = AseAtomsAdaptor.get_atoms(structures)
    view(atoms)
    return


