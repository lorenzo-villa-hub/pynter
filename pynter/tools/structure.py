#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 10:52:08 2020

@author: villa
"""

import numpy as np
import os.path as op
import os
import collections
from scipy.spatial import KDTree
from ase.visualize import view
from pymatgen.util.coord import pbc_shortest_vectors
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.core.trajectory import Trajectory

def _get_distance_vector_and_image(lattice,frac_coords1,frac_coords2,jimage=None):
    """
    Same as pymatgen.core.Lattice.get_distance_and_image but returns 
    the distance vetor instead of the norm.
    Gets distance between two frac_coords assuming periodic boundary
    conditions. If the index jimage is not specified it selects the j
    image nearest to the i atom and returns the distance and jimage
    indices in terms of lattice vector translations. If the index jimage
    is specified it returns the distance between the frac_coords1 and
    the specified jimage of frac_coords2, and the given jimage is also
    returned.
    Args:
        frac_coords1 (3x1 array): Reference fcoords to get distance from.
        frac_coords2 (3x1 array): fcoords to get distance from.
        jimage (3x1 array): Specific periodic image in terms of
            lattice translations, e.g., [1,0,0] implies to take periodic
            image that is one a-lattice vector away. If jimage is None,
            the image that is nearest to the site is found.
    Returns:
        (distance, jimage): distance and periodic lattice translations
        of the other site for which the distance applies. This means that
        the distance between frac_coords1 and (jimage + frac_coords2) is
        equal to distance.
    """
    if jimage is None:
        v, d2 = pbc_shortest_vectors(lattice, frac_coords1, frac_coords2, return_d2=True)
        fc = lattice.get_fractional_coords(v[0][0]) + frac_coords1 - frac_coords2  # type: ignore
        fc = np.array(np.round(fc), dtype=int)
        return v, fc

    jimage = np.array(jimage)
    mapped_vec = lattice.get_cartesian_coords(jimage + frac_coords2 - frac_coords1)  # type: ignore
    return mapped_vec, jimage  # type: ignore


def deform_lattice(structure,stdev=0.03):
    """
    Deform lattice vectors in structure. Also changes the angles of lattice vectors (off diagonal components).
    The function "apply_strain" from the Structure class only changes the base vectors' size but not the angles.

    Parameters
    ----------
    structure : 
        Structure object.
    stdev : float
        Standard deviation.

    Returns
    -------
    (Structure)
        Structure with deformed lattice.
    """
    def get_random_strain(stdev):
        D = np.eye(3) #strain matrix
        for i in range(3):
            D[i, i] += np.random.uniform(-stdev, stdev)
            for j in range(i + 1, 3):
                D[i, j] += np.random.uniform(-stdev, stdev) / 4
                D[j, i] = D[i, j]
        return D
    strain_matrix = get_random_strain(stdev=stdev)
    atoms = structure.to_ase_atoms()
    new_cell = atoms.get_cell().dot(strain_matrix)
    atoms.set_cell(new_cell,scale_atoms=True)
    return Structure.from_ase_atoms(atoms)


def get_distance_vector(site1,site2,jimage=None):
    """
    Get distance vector between two sites assuming periodic boundary conditions.
    Same as pymatgen.core.sites.Site.distance but returns a vector instead of norm.
    Args:
        other (PeriodicSite): Other site to get distance from.
        jimage (3x1 array): Specific periodic image in terms of lattice
            translations, e.g., [1,0,0] implies to take periodic image
            that is one a-lattice vector away. If jimage is None,
            the image that is nearest to the site is found.
    Returns:
        distance (float): Distance between the two sites
    """
    return _get_distance_vector_and_image(site1.lattice, site1.frac_coords, site2.frac_coords,jimage=jimage)[0]


def get_displacement_vectors(structures):
    """
    Get vectors in cartesian coords of site displacements w.r.t. the first structure.

    Parameters
    ----------
    structures : 
        list of Pymatgen Structure objects.

    Returns
    -------
    (numpy ndarray)
        Displacement vectors in cartesian coordinates.
    """
    traj = Trajectory.from_structures(structures,constant_lattice=True) #order of structures determines ref in traj
    traj.to_displacements()
    disp = traj.coords[1]
    return structures[0].lattice.get_cartesian_coords(disp)


def get_furthest_neighbors(structure, target_site, species):
    """
    Get the indexes of atomic sites that are furthest away from the target site
    and maximally spaced apart from each other, considering periodic boundary conditions.

    Parameters:
    - structure: pymatgen Structure object
    - target_site: a pymatgen PeriodicSite or index of the target site in the structure
    - species: list of species to find (ordered as given)

    Returns:
    - List of indices corresponding to the furthest sites of the given species
    """
    if isinstance(target_site, int):
        target_site = structure[target_site]

    selected_indices = []
    
    for specie in species:
        max_score = -1
        best_index = None

        for i, site in enumerate(structure):
            if site.species_string != specie or i in selected_indices:
                continue
            
            # Distance from target
            distance_to_target = structure.get_distance(structure.index(target_site),i)

            # Compute minimum distance to already selected sites
            if selected_indices:
                distances_to_selected = [structure.get_distance(j,i) for j in selected_indices]
                min_distance_to_selected = min(distances_to_selected)
            else:
                min_distance_to_selected = distance_to_target

            # Define score: prioritize sites far from target and far from other selections
            score = distance_to_target + min_distance_to_selected

            if score > max_score:
                max_score = score
                best_index = i
        
        if best_index is not None:
            selected_indices.append(best_index)

    return selected_indices

def get_WarrenCowley_order_parameter(
                        structure,
                        A_symbols,
                        X_symbols,
                        neighbors_cutoff):
    """
    Compute Warren Cowley order parameter to describe mixing tendencies in alloys

    The WC order parameter is computed as:
                alpha_{A-X} = 1 - < P_{A-X} > / C_X
    where:
    - <P_{A-X}> is the avg probability of finding X near A, computed as:
        <P_{A-X}> = n_X / n_total , n = number of nighbors of A
        C_X = concentration of X in structure
    
    Parameters
    ----------
    structure:
        Pymatgen Structure
    A_symbols: (list)
        List of element symbols to consider for A-site
    X_symbols: (list)
        List of element symbols to consider for X-site
    neighbors_cutoff: (float)
        Distance cutoff for neighbors counting

    Returns
    -------
    alpha: (float)
        Warren-Cowley order parameter
    """
    total_P_AX = 0
    n_A_sites = 0
    n_X_sites = 0
    for site in structure:
        if site.specie.symbol in A_symbols:
            n_A_sites += 1
            neighbors_list = structure.get_neighbors(site,r=neighbors_cutoff)
            n_X_neighbors = 0
            n_total_neighbors = 0
            for neighbor in neighbors_list:
                if neighbor.specie.symbol in X_symbols:
                    n_X_neighbors += 1
                    n_total_neighbors += 1
                elif neighbor.specie.symbol in A_symbols:
                    n_total_neighbors += 1
            total_P_AX += n_X_neighbors / n_total_neighbors
        elif site.specie.symbol in X_symbols:
            n_X_sites += 1

    if n_A_sites == 0 or n_X_sites == 0:
        raise ValueError("No A or X sites found to compute order parameter.")

    C_X = n_X_sites / (n_A_sites + n_X_sites)
    avg_P_AX = total_P_AX / n_A_sites
    alpha = 1 - avg_P_AX / C_X 
    return alpha


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
    is_site_in_structure,index = False,None
    check,index = is_site_in_structure_coords(site, structure,tol=tol)
    if check:
        s = structure[index]
        if site.specie.symbol == s.specie.symbol:
            is_site_in_structure =True
    
    return is_site_in_structure,index


def is_site_in_structure_coords(site, structure, tol=1e-3):
    """
    Check if a site's coordinates are present in a structure, using periodic boundary conditions.

    Parameters
    ----------
    site : (Site)
        PeriodicSite or Site object.
    structure : (Structure)
        Pymatgen Structure object.
    tol : (float), optional
        Tolerance for site comparison, normalized with respect to lattice size. 
        Default is 1e-03.

    Returns
    -------
    is_site_in_structure_coords : (bool)
        Whether the site exists in the structure within tolerance.
    index : (int or None)
        Index of matching site in structure if found, otherwise None.
    """
    # Normalize tolerance with lattice vector size
    l = site.lattice
    tol = np.sqrt(l.a**2 + l.b**2 + l.c**2) * tol

    # Convert structure sites to fractional coordinates
    frac_coords_structure = np.array([s.frac_coords for s in structure])
    kdtree_structure = KDTree(frac_coords_structure, boxsize=1.0)  # Take periodicity into account

    # Query the KDTree for nearest neighbor
    dist, index = kdtree_structure.query(site.frac_coords, distance_upper_bound=tol)

    if dist < tol:
        return True, index
    else:
        return False, None


def rattle_atoms(structure,stdev=0.05,seed=None):
    """
    Change randomly atomic positions. Uses the ASE.
    Similar to the method "perturb" from the Structure class
    Parameters
    ----------
    structure :
        Structure object
    stdev : float
        Standard deviation
    seed : int
        Seed for random number generation passed to the Atoms.rattle method.
    """
    atoms = structure.to_ase_atoms()
    atoms.rattle(stdev=stdev,seed=seed)
    return Structure.from_ase_atoms(atoms)


def remove_oxidation_state_from_site(site):
    """
    Remove oxidation state decoration from site.
    """
    new_sp = collections.defaultdict(float)
    for el, occu in site.species.items():
        sym = el.symbol
        new_sp[Element(sym)] += occu
    site.species = Composition(new_sp)
    return

def sort_sites_to_ref_coords(structure,structure_ref,extra_sites=[],tol=1e-03,get_indexes=False):
    """
    Sort Sites of one structure to match the order of coordinates in a reference structure. 

    Parameters
    ----------
    structure : (Structure)
        Structure to sort.
    structure_ref : (Structure)
        Reference Structure.
    extra_sites : (list), optional
        Sites to append at the end of the structure. The default is [].
    tol: (float)
        Tolerance for site comparison. The distance between sites in target and reference stucture is used, 
        periodicity is accounted for. The tolerance is normalized with respect to lattice vector size. 
        The default is 1e-03.
    get_indexes : (bool), optional
        Get list of mapping indexes for target structure sites in reference structure. The default is False.

    Returns
    -------
    new_structure : (Structure)
        Sorted Structure
    indexes : (list)
        If get_indexes is True a tuple is returned. List of mapping indexes for target structure sites
        in reference structure.
    """
    df = structure
    bk = structure_ref
    indexes = []
    new_sites=[]
    for s in df:
        check,index = is_site_in_structure_coords(s,bk,tol=tol)
        if check:
            indexes.append(index)
    for w in range(0,len(bk)):
        new_sites.insert(indexes[w],df[w])
    for s in extra_sites:
        new_sites.append(s)
            
    new_structure = df.copy()
    new_structure._sites = new_sites 
    if get_indexes:
        return new_structure, indexes
    else:
        return new_structure


def view_structures_with_ase(structures):
    """
    Visualize a Structure object or a list of Structure objects with the ASE. 
    First the Structure objects are converted into an ase Atom object, then "view" is used to visualize them.
    """
    if type(structures) == list:
        atoms=[]
        for s in structures:
            atoms.append(s.to_ase_atoms())
    else:
        atoms = structures.to_ase_atoms()
    view(atoms)
    return


def write_extxyz_file(file,structure,structure_ref=None,displacements=False):
    """
    Write extxyz format. Displacements can be included for visualization in OVITO. 

    Parameters
    ----------
    file : (str)
        Path to save extxyz file to.
    structure : (Structure)
        Structure to visualize.
    structure_ref : (Structure)
        Reference Structure if needed. The default is None.
    displacements : (bool)
        Include displacement vectors. The default is False.
    """
    atoms = structure.to_ase_atoms()
    if displacements:
        disp = get_displacement_vectors([structure_ref, structure])
        atoms.arrays.update({'disp':disp})
    if not op.exists(op.dirname(file)):
        os.makedirs(op.dirname(file))
    atoms.write(file,format='extxyz')
    return


def write_xdatcar_from_structures(structures,file='XDATCAR'):
    """
    Write XDATCAR file from a list of structures. The first structure determines the reference for the Trajectory object.
    """
    traj = Trajectory.from_structures(structures,constant_lattice=True) 
    if not op.exists(op.dirname(file)):
        os.makedirs(op.dirname(file))
    traj.write_Xdatcar(file)
    return
    
    
    
    
    
    
    
    
    
