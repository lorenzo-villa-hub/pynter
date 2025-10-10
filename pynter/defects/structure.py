#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 17:10:46 2023

@author: villa
"""

import warnings
import numpy as np
import os.path as op
import os

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
from pymatgen.core.trajectory import Trajectory

from pynter.tools.structure import is_site_in_structure_coords, sort_sites_to_ref_coords, write_extxyz_file
from .defects import Vacancy,Substitution,Interstitial,DefectComplex
from .generator import create_interstitials, create_vacancies, create_substitutions



"""Interstitial generator to be re-implemented using the new pymatgen defects"""
    
def create_interstitial_structures(structure,elements,supercell_size=None,**kwargs):
    """
    Create interstitial structures based on Voronoi with pymatgen,
    staring from a bulk structure (unit cell or supercell).
    Uses `Interstitial` objects generated with the `generator` module.

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
    structures : (list)
        List of interstitial structures
    """
    defects = create_interstitials(
                                structure=structure,
                                elements=elements,
                                supercell_size=supercell_size,
                                **kwargs)
    structures = [df.generate_defect_structure() for df in defects]
    return structures
    

def create_substitution_structures(structure,elements_to_replace,supercell_size=None,**kwargs):
    """
    Create substitution structures for each non-equivalent site,
    staring from a bulk structure (unit cell or supercell).
    Uses `Substitution` objects generated with the `generator` module.

    Parameters
    ----------
    structure : Structure
        Bulk structure, both unit cell or supercell can be used as input.
    elements_to_replace : (str)
        Dict with element symbol of specie to be replaced as keys and element 
        symbol of the species to be replaced with as values ({'old_El':'new_El'}).
    supercell_size : (int or numpy array)
        Input for the make_supercell function of the Structure class.
        If None the input structure is not modified. 
    kwargs : (dict)
        Kwargs to pass to SpaceGroupAnalyzer class.
        Args:
            symprec (float): Tolerance for symmetry finding. Defaults to 0.01,
            which is fairly strict and works well for properly refined
            structures with atoms in the proper symmetry coordinates. For
            structures with slight deviations from their proper atomic
            positions (e.g., structures relaxed with electronic structure
            codes), a looser tolerance of 0.1 (the value used in Materials
            Project) is often needed.
        angle_tolerance (float): Angle tolerance for symmetry finding. Defaults to 5 degrees.

    Returns
    -------
    structures : (list)
        List of substitution structures
    """
    defects = create_substitutions(
                                structure=structure,
                                elements_to_replace=elements_to_replace,
                                supercell_size=supercell_size,
                                **kwargs)
    structures = [df.generate_defect_structure() for df in defects]
    return structures


def create_vacancy_structures(structure,elements=None,supercell_size=None,**kwargs):
    """
    Create Vacancy objects for each non-equivalent site,
    staring from a bulk structure (unit cell or supercell).
    Uses `Vacancy` objects generated with the `generator` module.

    Parameters
    ----------
    structure : Structure
        Bulk structure, both unit cell or supercell can be used as input.
    elements : (str), optional
        Symbol of the elements for which vacancies are needed.
        If None all of the elements are considered. The default is None.
    supercell_size : (int or numpy array)
        Input for the make_supercell function of the Structure class.
        If None the input structure is not modified. 
    kwargs : (dict)
        Kwargs to pass to SpaceGroupAnalyzer class.
        Args:
            symprec (float): Tolerance for symmetry finding. Defaults to 0.01,
            which is fairly strict and works well for properly refined
            structures with atoms in the proper symmetry coordinates. For
            structures with slight deviations from their proper atomic
            positions (e.g., structures relaxed with electronic structure
            codes), a looser tolerance of 0.1 (the value used in Materials
            Project) is often needed.
        angle_tolerance (float): Angle tolerance for symmetry finding. Defaults to 5 degrees.

    Returns
    -------
    structures : (list)
        List of vacancy structures
    """
    defects = create_vacancies(
                            structure=structure,
                            elements=elements,
                            supercell_size=supercell_size,
                            **kwargs)
    structures = [df.generate_defect_structure() for df in defects]
    return structures

 
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
        if df_found.type=='DefectComplex':
            dfs = df_found.defects
        else:
            dfs = [df_found]
   
    for d in dfs:
        dsite = d.site
        dtype = d.type
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


def defect_finder(
                structure_defect,
                structure_bulk,
                tol=1e-3,
                max_number_of_defects=None,
                verbose=False):
    """
    Find defects by comparing defect and bulk structures.

    Parameters
    ----------
    structure_defect : Pymatgen Structure
        Defect structure.
    structure_bulk : Pymatgen Structure
        Bulk structure.
    tol : float
        Tolerance for fractional coordinates comparison (default is 1e-3).
    max_number_of_defects : (int)
        Impose a max number of defects to be found. If set, defects
        are ranked based on the coordinate distance btw defect and bulk sites
        (descending order), the first max_number_of_defects in the list 
        are given as the output (as single defect or defect complex).
    verbose : (bool)
        Print output.

    Returns
    -------
    Defect object (single or complex)
    """
    defects = []
    # Identify missing (vacancies) and additional (interstitials) sites
    matched_indices = set()
    for site in structure_defect:
        check, index = is_site_in_structure_coords(site,structure_bulk,tol=tol,return_distance=False)

        if check:
            matched_indices.add(index)  # Site exists in bulk, check for substitution
            if site.species != structure_bulk[index].species:
                defects.append(
                    Substitution(specie=site.specie.symbol,
                                defect_site=site,
                                bulk_structure=structure_bulk,
                                site_in_bulk=structure_bulk[index]))
        else:
            defects.append(Interstitial(specie=site.specie.symbol,
                                        defect_site=site,
                                        bulk_structure=structure_bulk))

    for j, site in enumerate(structure_bulk):
        if j not in matched_indices:
            defects.append(Vacancy(specie=site.specie.symbol,
                                        defect_site=site,
                                        bulk_structure=structure_bulk))

    if max_number_of_defects:
        filtered_defects = []
        defect_distances = []
        for df in defects:
            if df.type == 'Vacancy':
                check, i, d = is_site_in_structure_coords(df.site,structure_defect,tol=tol,return_distance=True)
                defect_distances.append((df,d))
            else:
                check, i, d = is_site_in_structure_coords(df.site,structure_bulk,tol=tol,return_distance=True)
                defect_distances.append((df,d))
        sorted_defect_distances = defect_distances.sort(key = lambda x : x[1], reverse=True)
        filtered_defects = [defect_distances[i][0] for i in range(max_number_of_defects)]
    else:
        filtered_defects = defects

    if len(filtered_defects) > 1:
        if len(filtered_defects) > 3:
            warnings.warn("More than 3 defects found, if not desired try to adjust the tolerance parameter ")
        defect = DefectComplex(filtered_defects, bulk_structure=structure_bulk)
    elif len(filtered_defects) == 1:
        defect = filtered_defects[0]
    else:
        warnings.warn("No defect has been found. Try to adjust the tolerance parameter.", UserWarning)
        defect = None
    
    if verbose:
        print(f'Defect automatically identified for defective structure with composition {structure_defect.composition}: \n {defect}')

    return defect


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