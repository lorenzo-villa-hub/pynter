#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:27:33 2020

@author: villa
"""

from pymatgen.analysis.defects.utils import StructureMotifInterstitial
import numpy as np
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun, Locpot, VolumetricData, Outcar
from pymatgen.analysis.defects.core import Vacancy , DefectEntry, Interstitial, Substitution, Defect
from pymatgen.analysis.defects.corrections import FreysoldtCorrection, KumagaiCorrection
import os.path as op
import importlib
from pynter.tools.structure import is_site_in_structure, is_site_in_structure_coords, sort_sites_to_ref_coords

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


def create_vacancy_structures(structure,elements=None,supercell_size=None):
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

    Returns
    -------
    vac_structures : (dict)
        Dictionary with vacancy types as keys and structures as values.

    """
    vac_structures={}
    bulk_structure = structure.copy()
    if supercell_size:
        bulk_structure.make_supercell(supercell_size)
    if not elements:
        elements = [el.symbol for el in bulk_structure.composition.elements]
    
    for el in bulk_structure.composition.elements:
        s = bulk_structure.copy()
        for site in s.sites:
            if el.symbol in elements:
                if site.specie == el:
                    s.remove_sites([bulk_structure.index(site)])
                    vac_structures[f'{el.symbol}'] = s
                    break

    return vac_structures


def create_substitution_structures(structure,replace,supercell_size=1):
    """
    Create structures with vacancies starting from a bulk structure (unit cell or supercell).

    Parameters
    ----------
    structure : Structure
        Bulk structure, both unit cell or supercell can be used as input.
    replace : (str), optional
        Dict with element symbol of specie to be replaced as keys and element 
        symbol of the species to be replaced with as values ({'old_El':{new_El}}).
    supercell_size : (int or numpy array), optional
        Input for the generate_defect_structure function of the Substitution class.

    Returns
    -------
    sub_structures : (dict)
        Dictionary with substitution types as keys and structures as values.

    """
    sub_structures={}
    bulk_structure = structure.copy()
    
    for el_to_sub in replace:
        for el in bulk_structure.composition.elements:
            s = bulk_structure.copy()
            for site in s.sites:
                if el.symbol == el_to_sub:
                    if site.specie == el:
                        sub_site = site
                        sub_el = replace[el.symbol]
                        defect_site = PeriodicSite(sub_el,sub_site.frac_coords,sub_site.lattice)
                        structure = Substitution(s,defect_site).generate_defect_structure(supercell_size)
                        sub_structures[f'{sub_el}-on-{el.symbol}'] = structure
                        break

    return sub_structures


def create_def_structure_for_visualization(structure_defect,structure_bulk,defects=None,sort_to_bulk=False,tol=1e-03):
    """
    Create defect structure for visualization in OVITO. The vacancies are shown by inserting 
    in the vacant site the element of same row and next group on the periodic table.
    If sort_to_bulk is True the Sites are sorted to match the Bulk structure.

    Parameters
    ----------
    structure_defect : (Pymatgen Structure object)
        Defect structure.
    structure_bulk : (Pymatgen Structure object)
        Bulk structure.
    defects : (tuple or list). 
        Tuple or list of tuples in the format (defect_site,defect_type)
        The format is the same as the output of defect_finder. If None defect_finder is used. The default is None.
    sort_to_bulk : (bool)
        Sort Sites of the defect structure to match the order of coordinates in the bulk structure
        (useful if the non-relaxed defect structure is not available). 
        If False only the dummy atoms are inserted and not further changes are made.
    tol: (float)
        Tolerance for site comparison. The distance between sites in target and reference stucture is used, 
        periodicity is accounted for. The tolerance is normalized with respect to lattice vector size. 
        The default is 1e-03.

    Returns
    -------
    new_structure (Pymatgen Structure object)
        Structure with dummy atoms as vacancies and interstitials in the bottom.
        The order of the Sites follow the order of the Bulk structure.

    """
    df = structure_defect.copy()
    bk = structure_bulk.copy()
    extra_sites=[]
    if defects:
        dfs = defects
    else:
        dfs = defect_finder(df,bk,tol=tol)
    
    # handle single defect case
    if not isinstance(dfs,list):
        dfs = [dfs] 
    for dsite,dtype in dfs:
        if dtype == 'Vacancy':
            check,i = is_site_in_structure_coords(dsite,bk,tol=tol)
            el = dsite.specie
            species = Element.from_row_and_group(el.row, el.group+1)
            df.insert(i=i,species=species,coords=dsite.frac_coords)
        elif dtype == 'Interstitial' and sort_to_bulk:
            extra_sites.append(dsite)
    # reorder to match bulk, useful if you don't have non-relaxed defect structure
    if sort_to_bulk:
        new_structure = sort_sites_to_ref_coords(df, bk, extra_sites,tol=tol)
    # In this case only dummy atoms are inserted, no further changes
    else:
        new_structure = df.copy()
    return new_structure
        


def defect_finder(structure_defect,structure_bulk,tol=1e-03):
    """
    Function to find defect comparing defect and bulk structure (Pymatgen objects). 

    Parameters
    ----------
    structure_defect : (Pymatgen Structure object)
        Defect structure.
    structure_bulk : (Pymatgen Structure object)
        Bulk structure.
    tol: (float)
        Tolerance for fractional coordinates comparison.

    Returns
    -------
    defects : (list)
        List of tuples with defect site (Pymatgen PeriodicSite object) and defect type ("Vacancy","Interstitial" or "Substitution")
        for each point defect present. If len(defects) is 1 only the tuple is returned.
    """
    df = structure_defect
    bk = structure_bulk
    defects = []
    
    for s in df:
        if is_site_in_structure_coords(s,bk,tol=tol)[0]:
            if is_site_in_structure(s,bk,tol=tol)[0] == False:
                dsite = s
                dtype = 'Substitution'
                defects.append((dsite,dtype))
        else:
            dsite = s
            dtype = 'Interstitial'
            defects.append((dsite,dtype))
            
    for s in bk:
        if is_site_in_structure_coords(s,df,tol=tol)[0] == False:
            dsite = s
            dtype = 'Vacancy'
            defects.append((dsite,dtype))
    
    if len(defects) == 1:
        return defects[0]
    else:
        return defects
            
    
    
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
    
    structure_bulk = locpot_pure.structure
    defect_site = PeriodicSite(defect_specie, coords=defect_site_coordinates, lattice = locpot_pure.structure.lattice)
    
    module = importlib.import_module("pymatgen.analysis.defects.core")
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
        defect_site, defect_type = defect_finder(structure_defect, structure_bulk, tol=tol)
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
    
    module = importlib.import_module("pymatgen.analysis.defects.core")
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
    
    
def get_kumagai_correction_from_jobs(job_defect,job_bulk,dielectric_tensor,defect_site=None,sampling_radius=None,
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
    defect_site : (Site) , optional
        Defect site. If None the defect site is found unsing defect_finder. The default is None.
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
    
    if not defect_site:
        defect_site, defect_type = defect_finder(structure_defect, structure_bulk,tol)
        defect_specie = defect_site.specie.symbol
    
    charge = job_defect.charge
    
    corr = get_kumagai_correction(structure_defect, structure_bulk, path_to_defect_outcar, path_to_bulk_outcar,
                                  dielectric_tensor, charge, defect_type, defect_specie, defect_site, sampling_radius,
                                  gamma, tol, get_plot)
    
    return corr
        
    
    
    