#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:35:23 2021

@author: villa
"""

import numpy as np
import pandas as pd
import matplotlib



def get_elastic_dipole_tensor(stress_defect,stress_bulk,bulk_volume,corrections={}):
    """
    Calculate elastic dipole tensor from stresses.

    Parameters
    ----------
    stress_defect : (np.array)
        Stresses from defect calculation in kbar (VASP output)
    stress_bulk : (np.array)
        Stresses from bulk calculation in kbar (VASP output)
    bulk_volume : (float)
        Cell volume of bulk calculation in A°^3.
    corrections : (bool)
        Add correction terms to the residual stress tensor.

    Returns
    -------
    dipole_tensor : (np.array)
        Elastic dipole tensor in eV.
    """    
    res_stress = np.array(stress_defect) - np.array(stress_bulk)
    dipole_tensor = -1*bulk_volume*res_stress #sign is inverted with respect to VASP output
    dipole_tensor += sum([v for k,v in corrections.items()]) 
    return dipole_tensor


def get_relaxation_volume(stress_defect,stress_bulk,bulk_modulus,bulk_volume,corrections={}):
    """
    Calculate relaxation volume from stresses.

    Parameters
    ----------
    stress_defect : (np.array)
        Stresses from defect calculation in kbar (VASP output)
    stress_bulk : (np.array)
        Stresses from bulk calculation in kbar (VASP output)
    bulk_volume : (float)
        Cell volume of bulk calculation in A°^3.
    bulk_modulus : (float)
        Bulk modulus in GPa.
    corrections : (bool)
        Add correction terms to the residual stress tensor.
        
    Returns
    -------
    rel_volume : (float)
        Relaxation volume in A°^3.
    """
    bulk_modulus = bulk_modulus*10 # from GPa to kbar
    dipole_tensor = get_elastic_dipole_tensor(
                                            stress_defect=stress_defect,
                                            stress_bulk=stress_bulk,
                                            bulk_volume=bulk_volume,
                                            corrections=corrections)
    
    pressure = np.trace(dipole_tensor)/3
    rel_volume = pressure/bulk_modulus
    return rel_volume
