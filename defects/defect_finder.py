#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 19:08:20 2020

@author: lorenzo
"""

def defect_finder(structure_defect,structure_bulk):
    """
    Function to find defect comparing defect and bulk structure (Pymatgen objects).
    Warning: apparantely comparing a structure read from Vasprun and one read from a Poscar doesn't work. 

    Parameters
    ----------
    structure_defect : (Pymatgen Structure object)
        Defect structure.
    structure_bulk : (Pymatgen Structure object)
        Bulk structure.

    Returns
    -------
    defect_site : (Pymatgen PeriodicSite object)
    defect_type : (str)
        Type of defect ("Vacancy","Interstitial" or "Substitution").
    """

    df = structure_defect
    bk = structure_bulk
    
    if len(df) > len(bk):
        defect_type = 'Interstitial'
        for s in df:
            if s not in bk:
                defect_site = s
        
    elif len(df) < len(bk):
        defect_type = 'Vacancy'
        for s in bk:
            if s not in df:
                defect_site = s

    elif len(df) == len(bk):
        comp_df = structure_defect.composition
        comp_bk = structure_bulk.composition
        if comp_df != comp_bk:
            defect_type = 'Substitution'
            for el in comp_df:
                if el not in comp_bk:
                    sub_element = el
            for s in df:
                if s.specie == sub_element:
                    defect_site = s

    return defect_site  , defect_type                  