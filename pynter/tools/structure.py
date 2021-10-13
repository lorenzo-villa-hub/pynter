#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 10:52:08 2020

@author: villa
"""

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view


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
    is_site_in_structure = False
    for s in structure:
        if np.allclose(site.frac_coords,s.frac_coords,atol=tol) and site.specie.symbol == s.specie.symbol:
            is_site_in_structure = True
            index = structure.index(s)
            return is_site_in_structure,index
    index=None
    return is_site_in_structure,index


def is_site_in_structure_coords(site,structure,tol=1e-03):
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


