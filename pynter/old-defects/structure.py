#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 17:10:46 2023

@author: villa
"""

import warnings
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
import os.path as op
import os
from pynter.tools.structure import is_site_in_structure_coords, sort_sites_to_ref_coords, write_extxyz_file
from pynter.defects.defects import Vacancy,Substitution,Interstitial,DefectComplex
from pymatgen.core.trajectory import Trajectory

"""Interstitial generator to be re-implemented using the new pymatgen defects"""
    
def create_interstitial_structures(structure,elements,supercell_size=None,**kwargs):
    """
    Create structures with interstitials based on Voronoi with pymatgen. 

    Parameters
    ----------
    structure : (Structure)
        Bulk structure.
    elements : (list)
        List of element symbols.
    supercell_size : (int), optional
        Input for the make_supercell function of the Structure class.
        If None the input structure is not modified. The default is None.
    kwargs: 
        Arguments to pass to VoronoiInterstitialGenerator:
            clustering_tol: Tolerance for clustering the Voronoi nodes.
            min_dist: Minimum distance between an interstitial and the nearest atom.
            ltol: Tolerance for lattice matching.
            stol: Tolerance for structure matching.
            angle_tol: Angle tolerance for structure matching.
            kwargs: Additional keyword arguments for the ``TopographyAnalyzer`` constructor.

    Returns
    -------
    interstitial_structures : (dict)
        Dictionary with element as keys and list of structures as values.
    """
    from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
    bulk_structure = structure.copy()
    if supercell_size:
        bulk_structure.make_supercell(supercell_size)
    generator = VoronoiInterstitialGenerator().generate(bulk_structure,elements)
    interstitial_structures = {}
    for inter in generator:
        el = inter.site.specie.element.symbol
        if el not in interstitial_structures.keys():
            interstitial_structures[el] = []
        interstitial_structures[el].append(inter.defect_structure)
            
    return interstitial_structures
    

def create_substitution_structures(structure,replace,supercell_size=1):
    """
    Create structures with substitutions starting from a bulk structure (unit cell or supercell).

    Parameters
    ----------
    structure : Structure
        Bulk structure, both unit cell or supercell can be used as input.
    replace : (str), optional
        Dict with element symbol of specie to be replaced as keys and element 
        symbol of the species to be replaced with as values ({'old_El':{new_El}}).
    supercell_size : (int or numpy array), optional
        Input for the generate_defect_structure function of the old pymatgen Substitution class.

    Returns
    -------
    sub_structures : (dict)
        Dictionary with substitution types as keys and structures as values.
    """
    from pynter.defects.pmg.pmg_defects_core import Substitution
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
        List of defect objects.
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
        df_found = defect_finder(df,bk,tol=tol)
        if df_found.defect_type=='DefectComplex':
            dfs = df_found.defects
        else:
            dfs = [df_found]
   
    for d in dfs:
        dsite = d.site
        dtype = d.defect_type
        if dtype == 'Vacancy':
            check,i = is_site_in_structure_coords(dsite,bk,tol=tol)
            el = dsite.specie
            species = Element.from_row_and_group(el.row, el.group+1)
            df.insert(idx=i,species=species,coords=dsite.frac_coords)
        elif dtype == 'Interstitial' and sort_to_bulk:
            extra_sites.append(dsite)

    # reorder to match bulk, useful if you don't have the non-relaxed defect structure
    if sort_to_bulk:
        new_structure = sort_sites_to_ref_coords(df, bk, extra_sites,tol=tol)
    # In this case only dummy atoms are inserted, no further changes
    else:
        new_structure = df.copy()
    return new_structure        


def defect_finder(structure_defect, structure_bulk, tol=1e-3):
    """
    Optimized function to find defects by comparing defect and bulk structures.

    Parameters
    ----------
    structure_defect : Pymatgen Structure
        Defect structure.
    structure_bulk : Pymatgen Structure
        Bulk structure.
    tol : float, optional
        Tolerance for fractional coordinates comparison (default is 1e-3).

    Returns
    -------
    Defect object (single or complex)
    """
    defects = []
    # Identify missing (vacancies) and additional (interstitials) sites
    matched_indices = set()
    for site in structure_defect:
        check, index = is_site_in_structure_coords(site,structure_bulk,tol=tol)
        
        if check:
            matched_indices.add(index)  # Site exists in bulk, check for substitution
            if site.species != structure_bulk[index].species:
                defects.append(
                    Substitution(site,structure_bulk,site_in_bulk=structure_bulk[index])
                    )
        else:
            defects.append(Interstitial(site, structure_bulk))

    for j, site in enumerate(structure_bulk):
        if j not in matched_indices:
            defects.append(Vacancy(site, structure_bulk))

    if len(defects) > 1:
        if len(defects) > 3:
            warnings.warn("More than 3 defects found, if not desired try to adjust the tolerance parameter ")
        return DefectComplex(defects, structure_bulk)
    elif len(defects) == 1:
        return defects[0]
    else:
        warnings.warn("No defect has been found. Try to adjust the tolerance parameter.", UserWarning)
        return defects



def get_trajectory_for_visualization(structure_defect,structure_bulk,defects=None,tol=1e-03,file=None):
    """
    Create trajectory from defect and bulk structures for visualization in OVITO. 
    The vacancies are shown by inserting in the vacant site the element of same row and next group on the periodic table.

    Parameters
    ----------
    structure_defect : (Pymatgen Structure object)
        Defect structure.
    structure_bulk : (Pymatgen Structure object)
        Bulk structure.
    defects : (tuple or list). 
        Tuple or list of tuples in the format (defect_site,defect_type)
        The format is the same as the output of defect_finder. If None defect_finder is used. The default is None.
    tol: (float)
        Tolerance for site comparison. The distance between sites in target and reference stucture is used, 
        periodicity is accounted for. The tolerance is normalized with respect to lattice vector size. 
        The default is 1e-03.
    file : (str)
        File to save XDATCAR. 
    Returns
    -------
    new_structure (Pymatgen Structure object)
        Structure with dummy atoms as vacancies and interstitials in the bottom.
        The order of the Sites follow the order of the Bulk structure.

    """
    sb = structure_bulk
    dummy = create_def_structure_for_visualization(structure_defect, structure_bulk,defects,sort_to_bulk=True,tol=tol)
    traj = Trajectory.from_structures([dummy,sb],constant_lattice=True)     
    if file:
        if not op.exists(op.dirname(file)):
            os.makedirs(op.dirname(file))
        traj.write_Xdatcar(file)
    return traj
        
        
def write_extxyz_for_visualization(file,structure_defect,structure_bulk,defects=None,tol=1e-03): 
    """
    Write extxyz file for visualization in OVITO. The displacements w.r.t the bulk structure are included.
    The vacancies are shown by inserting in the vacant site the element of same row and next group on the periodic table.
    
    Parameters
    ----------
    file : (str)
        Path to save file. 
    structure_defect : (Pymatgen Structure object)
        Defect structure.
    structure_bulk : (Pymatgen Structure object)
        Bulk structure.
    defects : (tuple or list). 
        Tuple or list of tuples in the format (defect_site,defect_type)
        The format is the same as the output of defect_finder. If None defect_finder is used. The default is None.
    tol: (float)
        Tolerance for site comparison. The distance between sites in target and reference stucture is used, 
        periodicity is accounted for. The tolerance is normalized with respect to lattice vector size. 
        The default is 1e-03.
    """
    sb = structure_bulk
    dummy = create_def_structure_for_visualization(structure_defect, structure_bulk,defects,sort_to_bulk=True,tol=tol)
    write_extxyz_file(file, dummy, sb,displacements=True)
    return