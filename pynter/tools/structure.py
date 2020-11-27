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
        if np.allclose(site.frac_coords,s.frac_coords,rtol=tol) and site.specie.symbol == s.specie.symbol:
            is_site_in_structure = True
            index = structure.index(s)
            return is_site_in_structure,index
    index=None
    return is_site_in_structure,index


def view_structure_with_ase(structure):
    """
    Visualize the Structure object with the ASE. First the Structure object is converted into an
    ase Atom object, the "view" is used to visualize it.
    """
    atoms = AseAtomsAdaptor.get_atoms(structure)
    view(atoms)
    return