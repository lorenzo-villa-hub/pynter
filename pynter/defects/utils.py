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
from pymatgen.io.vasp.outputs import Vasprun, Locpot, VolumetricData, Outcar
from pymatgen.analysis.defects.core import Vacancy , DefectEntry, Interstitial, Substitution, Defect
from pymatgen.analysis.defects.corrections import FreysoldtCorrection, KumagaiCorrection
import os.path as op
import importlib

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


def get_delta_atoms(structure_defect,structure_bulk):
    """
    Function to biuld delta_atoms dictionary starting from Pymatgen Structure objects.
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
    delta_atoms = {}
    comp_defect = structure_defect.composition
    comp_bulk = structure_bulk.composition
    for el,n in comp_defect.items():
        nsites_defect = n
        nsites_bulk = comp_bulk[el] if el in comp_bulk.keys() else 0
        delta_n = nsites_defect - nsites_bulk
        if delta_n != 0:
            delta_atoms[el] = delta_n
        
    return delta_atoms



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



def get_kumagai_correction(structure_defect,structure_bulk,path_to_defect_outcar,path_to_bulk_outcar,dielectric_tensor,
                           charge,defect_type=None,defect_specie=None,defect_site=None,sampling_radius=None,gamma=None,
                           get_plot=False):
    
    
    if not defect_site and not defect_type and not defect_specie:
        defect_site, defect_type = defect_finder(structure_defect, structure_bulk)
        defect_specie = defect_site.specie.symbol
   
    site_matching_indices = []
    for site in structure_defect:
        site_matching_indices.append([structure_bulk.index(site),structure_defect.index(site)])
    # for sb in structure_bulk:
    #     if sb != defect_site:
    #         for sd in structure_defect:
    #             if sb == sd:
    #                 site_matching_indices.append([structure_bulk.index(sb),structure_bulk.index(sd)])
    
    print(site_matching_indices)
    bulk_atomic_site_averages = Outcar(op.join(path_to_bulk_outcar,'OUTCAR')).read_avg_core_poten()[-1]
    defect_atomic_site_averages = Outcar(op.join(path_to_defect_outcar,'OUTCAR')).read_avg_core_poten()[0]
    defect_frac_sc_coords = defect_site.frac_coords
    initial_defect_structure = structure_defect
    
    parameters = {}
    parameters['bulk_atomic_site_averages'] = bulk_atomic_site_averages
    parameters['defect_atomic_site_averages'] = defect_atomic_site_averages
    parameters['site_matching_indices'] = site_matching_indices
    parameters['initial_defect_structure'] = initial_defect_structure
    parameters['defect_frac_sc_coords'] = defect_frac_sc_coords
    
    module = importlib.import_module("pymatgen.analysis.defects.core")
    defect_class = getattr(module,defect_type)
    defect = defect_class(structure_bulk, defect_site, charge=charge, multiplicity=None)
    defect_entry = DefectEntry(defect,None,corrections=None,parameters=parameters)

    kumagai = KumagaiCorrection(dielectric_tensor,sampling_radius,gamma)
    kumagai_corrections = kumagai.get_correction(defect_entry)
    
    if get_plot:
        plt = kumagai.plot(1)
        return kumagai_corrections , plt
    else:    
        return kumagai_corrections
    
    
def get_kumagai_correction_from_jobs(job_defect,job_bulk,dielectric_tensor,defect_site=None,sampling_radius=None,
                                     gamma=None,get_plot=False):

    structure_defect = job_defect.initial_structure
    structure_bulk = job_bulk.final_structure
    path_to_defect_outcar = op.join(job_defect.path)
    path_to_bulk_outcar = op.join(job_bulk.path)
    
    if not defect_site:
        defect_site, defect_type = defect_finder(structure_defect, structure_bulk)
        defect_specie = defect_site.specie.symbol
    
    charge = job_defect.charge
    
    corr = get_kumagai_correction(structure_defect, structure_bulk, path_to_defect_outcar, path_to_bulk_outcar,
                                  dielectric_tensor, charge, defect_type, defect_specie, defect_site, sampling_radius,
                                  gamma, get_plot)
    
    return corr
    
    
    
    
    
    
    