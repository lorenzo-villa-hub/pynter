#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:27:33 2020

@author: villa
"""

from pymatgen.core.periodic_table import Element

def convert_conc_from_weight_to_cm3(c_weight,target_el,composition,bulk_structure):
    """
    Convert concentration of dopand from weight % (common experimental data) to cm^-3.

    Parameters
    ----------
    c_weight : (float)
        Defect concentration in weight %.
    target_el : (str)
        Symbol of target element.
    composition : (Composition)
        Composition (Pymatgen) of the material, uses only list of elements.
    bulk_structure : (Structure)
        Structure (Pymatgen) of bulk, needed for lattice volume.

    Returns
    -------
    conc : (float)
        Concentration in cm-3^.

    """
    sum_MM = sum([el.atomic_mass for el in composition.elements])
    r_MM = Element(target_el).atomic_mass / sum_MM
    conc = c_weight/r_MM * 1/bulk_structure.lattice.volume * 1e24 *0.01
    return conc




     