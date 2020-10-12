#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:27:33 2020

@author: villa
"""

from pymatgen.analysis.defects.utils import StructureMotifInterstitial
import numpy as np
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun, Locpot, VolumetricData
from pymatgen.analysis.defects.core import Vacancy , DefectEntry, Interstitial, Substitution
from pymatgen.analysis.defects.corrections import FreysoldtCorrection

def create_interstitial_supercells(structure,element,size=2):
    
    '''
    Create interstitial structures with InFiT (Interstitialcy Finding Tool) algorithm with Pymatgen.
    
    Parameters
    ----------
    structure: (Pymatgen Structure)
        Bulk structure.
    element: (str or Element)
        interstitial element.
    
    Returns
    -------
    interstitials: (dict)
        Dictionary with sites (Site object) and supercell structures (Structure object)
    '''
    
    int_object = StructureMotifInterstitial(structure,element)
    int_sites = int_object.enumerate_defectsites()
    int_supercells = int_object.make_supercells_with_defects([[size,0,0],[0,size,0],[0,0,size]])
    
    interstitials = {}
    interstitial_sites =  {}
    interstitial_structures = {}

    for i in range(0,len(int_sites)):        
        int_coord_type = int_object.get_motif_type(i)
        interstitial_sites[int_coord_type] = int_sites[i]
        struct_int = int_supercells[i+1]
        interstitial_structures[int_coord_type] = struct_int
        
        interstitials['sites'] = interstitial_sites
        interstitials['structures'] = interstitial_structures
        
    return interstitials


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



def get_freysoldt_correction(defect_type, defect_specie, path_to_defect_locpot,path_to_pure_locpot,charge,
                             dielectric_constant,defect_site_coordinates,energy_cutoff=500,get_plot=False):
    
    ''' Function to perform charge corrections according to the method proposed py Freysoldt
        If this correction is used, please reference Freysoldt's original paper.
        doi: 10.1103/PhysRevLett.102.016402
        
        Args:
            defect_type: 'vacancy' or 'interstitial'
            defect_specie: string with element occupying the defect site
            path_to_defect_locpot: path to LOCPOT file of defect structure
            path_to_pure_locpot: path to LOCPOT file of Pure structure
            charge: Charge of the defected system
            dielectric_constant: Dielectric constant
            defect_site_coordinates: numpy array with fractional coordinates of defect site
            energy_cutoff: Cut-off of plane wave expansion
            get_plot: return also Matplotlib object with plot
            
        Returns:
            Freysoldt corrections values as a dictionary 
            '''
    # acquiring data from LOCPOT files    
    locpot_pure = Locpot.from_file(path_to_pure_locpot)
    vol_data_pure = VolumetricData(locpot_pure.structure,locpot_pure.data)
    
    locpot_defect = Locpot.from_file(path_to_defect_locpot)
    vol_data_defect = VolumetricData(locpot_defect.structure,locpot_defect.data)
    
    parameters = {}
    parameters['axis_grid'] = []
    parameters['bulk_planar_averages'] = []
    parameters['defect_planar_averages'] = []
    for i in range(0,3):
        parameters['axis_grid'].append(vol_data_pure.get_axis_grid(i))
        parameters['bulk_planar_averages'].append(vol_data_pure.get_average_along_axis(i))
        parameters['defect_planar_averages'].append(vol_data_defect.get_average_along_axis(i))
    parameters['initial_defect_structure'] = locpot_defect.structure
    parameters['defect_frac_sc_coords'] = defect_site_coordinates
    
    defect_site = PeriodicSite(defect_specie, coords=defect_site_coordinates, lattice = locpot_pure.structure.lattice)
    
    if defect_type == 'Vacancy':
        defect = Vacancy(locpot_pure.structure, defect_site, charge=charge, multiplicity=None)
    if defect_type == 'Interstitial':
        defect = Interstitial(locpot_pure.structure, defect_site, charge=charge, multiplicity=None)
    if defect_type == 'Substitution':
        defect = Substitution(locpot_pure.structure, defect_site, charge=charge, multiplicity=None)
        
    defect_entry = DefectEntry(defect,None,corrections=None,parameters=parameters)
    
    freysoldt_class = FreysoldtCorrection(dielectric_constant,energy_cutoff=energy_cutoff)
    
    freysoldt_corrections = freysoldt_class.get_correction(defect_entry)
  
    if get_plot:
        plt = freysoldt_class.plot(1)
        return freysoldt_corrections , plt
    else:    
        return freysoldt_corrections







