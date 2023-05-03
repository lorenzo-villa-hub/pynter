#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:26:50 2023

@author: villa
"""
from pymatgen.io.vasp.outputs import Vasprun, Locpot, VolumetricData, Outcar
from pynter.defects.pmg.pmg_defects_core import Vacancy , DefectEntry, Interstitial, Substitution
from pynter.defects.pmg.pmg_defects_corrections import FreysoldtCorrection, KumagaiCorrection
import os.path as op
import importlib
from pymatgen.core.sites import PeriodicSite
from pynter.tools.structure import is_site_in_structure_coords
from pynter.defects.structure import defect_finder


def get_freysoldt_correction(defect_type, defect_specie, path_to_defect_locpot,path_to_pure_locpot,charge,
                             dielectric_constant,defect_site_coordinates,energy_cutoff=500,get_plot=False):
    
    ''' Function to perform charge corrections according to the method proposed py Freysoldt
        If this correction is used, please reference Freysoldt's original paper.
        doi: 10.1103/PhysRevLett.102.016402
        
        Args:
            defect_type: 'Vacancy', 'Interstitial' or 'Substitution'
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
    
    structure_bulk = locpot_pure.structure
    defect_site = PeriodicSite(defect_specie, coords=defect_site_coordinates, lattice = locpot_pure.structure.lattice)
    
    module = importlib.import_module("pynter.defects.pmg.pmg_defects_core")
    defect_class = getattr(module,defect_type)
    defect = defect_class(structure_bulk, defect_site, charge=charge, multiplicity=None)
    defect_entry = DefectEntry(defect,None,corrections=None,parameters=parameters)
    
    freysoldt_class = FreysoldtCorrection(dielectric_constant,energy_cutoff=energy_cutoff)
    
    freysoldt_corrections = freysoldt_class.get_correction(defect_entry)
  
    if get_plot:
        plt = freysoldt_class.plot(1)
        return freysoldt_corrections , plt
    else:    
        return freysoldt_corrections


def get_freysoldt_correction_from_jobs(job_defect,job_bulk,dielectric_constant,energy_cutoff=500,tol=1e-03,get_plot=False):
    """
    Get Freysoldt corrections from VaspJob objects.

    Parameters
    ----------
    job_defect : (VaspJob)
        Defect calculation.
    job_bulk : (VaspJob)
        Bulk calculation.
    dielectric_constant : (float)
        Dielectric constant
    energy_cutoff : (int)
        Cut-off of plane wave expansion in eV
    tol : (float)
        Tolerance for defect_finder function. The default is 1e-03.
    get_plot : (bool), optional
        Get Matplotlib object with plot. The default is False.

    Returns
    -------
    corr : (dict or tuple)
        Dictionary with corrections, if get_plot is True a tuple with dict and plt object is returned.

    """
    structure_defect = job_defect.initial_structure
    structure_bulk = job_bulk.final_structure
    path_to_defect_locpot = op.join(job_defect.path,'LOCPOT')
    path_to_pure_locpot = op.join(job_bulk.path,'LOCPOT')

    df_found = defect_finder(structure_defect, structure_bulk, tol=tol)
    defect_site, defect_type = df_found.site, df_found.defect_type
    defect_specie = defect_site.specie.symbol
    defect_site_coordinates = defect_site.frac_coords
    charge = job_defect.charge

    corr = get_freysoldt_correction(defect_type, defect_specie, path_to_defect_locpot, path_to_pure_locpot, 
                                    charge, dielectric_constant, defect_site_coordinates,energy_cutoff,get_plot)
    
    return corr



def get_kumagai_correction(structure_defect,structure_bulk,path_to_defect_outcar,path_to_bulk_outcar,dielectric_tensor,
                           charge,defect_type=None,defect_specie=None,defect_site=None,sampling_radius=None,gamma=None,
                           tol=1e-03,get_plot=False):
    """
    Get Kumagai correction with Pymatgen.

    Parameters
    ----------
    structure_defect : (Structure)
        Structure of defect.
    structure_bulk : (Structure)
        Bulk structure.
    path_to_defect_outcar : (str)
        Path to OUTCAR of defect calculation.
    path_to_bulk_outcar : (str)
        Path to OUTCAR of pure calculation.
    dielectric_tensor : (array or float)
        Dielectric tensor, if is a float a diagonal matrix is constructed.
    charge : (int or float)
        Charge of the defect.
    defect_type : (str), optional
        Type of defect ('Vacancy','Interstitial' or 'Substitution')
        If None it's determined with defect_finder. The default is None.
    defect_specie : (str), optional
        Symbol of the defect specie.
        If None it's determined with defect_finder. The default is None.
    defect_site : (Site), optional
        Site of defect. If None it's determined with defect_finder. The default is None.
    sampling_radius (float): radius (in Angstrom) which sites must be outside
        of to be included in the correction. Publication by Kumagai advises to
        use Wigner-Seitz radius of defect supercell, so this is default value.
    gamma (float): convergence parameter for gamma function.
                    Code will automatically determine this if set to None.
    tol : (float)
        Tolerance for comparing sites and defect_finder function. The default is 1e-03.
    get_plot : (bool), optional
        Get Matplotlib object with plot. The default is False.

    Returns
    -------
    corr : (dict or tuple)
        Dictionary with corrections, if get_plot is True a tuple with dict and plt object is returned.
    """
    
    if not defect_site and not defect_type and not defect_specie:
        df_found = defect_finder(structure_defect, structure_bulk, tol=tol)
        defect_site, defect_type = df_found.site, df_found.defect_type
        defect_specie = defect_site.specie.symbol
   
    site_matching_indices = []
    for site in structure_defect:
        site_in_str ,index_bulk = is_site_in_structure_coords(site, structure_bulk,tol=tol)
        if site_in_str:
            site_matching_indices.append([index_bulk,structure_defect.index(site)])
        else:
            print(f'Warning in Kumagai corrections: Site {site} is not in bulk structure')
    
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
    
    module = importlib.import_module("pynter.defects.pmg.pmg_defects_core")
    defect_class = getattr(module,defect_type)
    defect = defect_class(structure_bulk, defect_site, charge=charge, multiplicity=1)
    defect_entry = DefectEntry(defect,None,corrections=None,parameters=parameters)

    kumagai = KumagaiCorrection(dielectric_tensor,sampling_radius,gamma)
    kumagai_corrections = kumagai.get_correction(defect_entry)
    
    if get_plot:
        plt = kumagai.plot()
        return kumagai_corrections , plt
    else:    
        return kumagai_corrections
    
    
def get_kumagai_correction_from_jobs(job_defect,job_bulk,dielectric_tensor,sampling_radius=None,
                                     gamma=None,tol=1e-03,get_plot=False):
    """
    Get Kumagai corrections from VaspJob objects.

    Parameters
    ----------
    job_defect : (VaspJob)
        Defect calculation.
    job_bulk : (VaspJob)
        Bulk calculation.
    dielectric_tensor : (array or float)
        Dielectric tensor, if is a float a diagonal matrix is constructed.
    sampling_radius (float): radius (in Angstrom) which sites must be outside
        of to be included in the correction. Publication by Kumagai advises to
        use Wigner-Seitz radius of defect supercell, so this is default value.
    gamma (float): convergence parameter for gamma function.
                    Code will automatically determine this if set to None.
    tol : (float)
        Tolerance for defect_finder function. The default is 1e-03.
    get_plot : (bool), optional
        Get Matplotlib object with plot. The default is False.

    Returns
    -------
    corr : (dict or tuple)
        Dictionary with corrections, if get_plot is True a tuple with dict and plt object is returned.

    """
    
    structure_defect = job_defect.initial_structure
    structure_bulk = job_bulk.final_structure
    path_to_defect_outcar = op.join(job_defect.path)
    path_to_bulk_outcar = op.join(job_bulk.path)

    df_found = defect_finder(structure_defect, structure_bulk, tol=tol)
    defect_site, defect_type = df_found.site, df_found.defect_type
    defect_specie = defect_site.specie.symbol
    
    charge = job_defect.charge
    
    corr = get_kumagai_correction(structure_defect, structure_bulk, path_to_defect_outcar, path_to_bulk_outcar,
                                  dielectric_tensor, charge, defect_type, defect_specie, defect_site, sampling_radius,
                                  gamma, tol, get_plot)
    
    return corr
