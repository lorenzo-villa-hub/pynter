
from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pynter.defects.defects import Interstitial, Substitution, Vacancy
from pynter.tools.structure import remove_oxidation_state_from_site



def create_interstitials(structure,elements,supercell_size=None,**kwargs):
    """
    Create Interstitial objects based on Voronoi with pymatgen,
    staring from a bulk structure (unit cell or supercell).

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
    defects : (list)
        List of Interstitial objects
    """
    from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
    
    defects = []
    bulk_structure = structure.copy()
    if supercell_size:
        bulk_structure.make_supercell(supercell_size)
    generator = VoronoiInterstitialGenerator().generate(bulk_structure,elements)
    for inter in generator:
        bulk_structure.remove_oxidation_states()
        remove_oxidation_state_from_site(inter.site)
        interstitial = Interstitial(
                                defect_site=inter.site,
                                bulk_structure=bulk_structure,
                                multiplicity=inter.multiplicity,
                                label=f'mult{inter.multiplicity}')
        defects.append(interstitial)
    return defects


def create_substitutions(structure,elements_to_replace,supercell_size=None, **kwargs):
    """
    Create Substitution objects for each non-equivalent site,
    staring from a bulk structure (unit cell or supercell).

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
    defects : (list)
        List of Substitution objects
    """
    defects = [] 
    bulk_structure = structure.copy()
    if supercell_size:
        bulk_structure.make_supercell(supercell_size)
    else:
        sites = bulk_structure.sites

    sym_struct = SpacegroupAnalyzer(structure=bulk_structure,**kwargs).get_symmetrized_structure()
    for el_to_sub,el_subbed in elements_to_replace.items():
        for site_group in sym_struct.equivalent_sites:
            site = site_group[0]
            if site.specie.symbol == el_to_sub:
                defect_site = PeriodicSite(el_subbed,site.frac_coords,site.lattice)
                substitution = Substitution(
                                        defect_site=defect_site,
                                        bulk_structure=bulk_structure,
                                        site_in_bulk=site)
                multiplicity = len(sym_struct.find_equivalent_sites(site))
                substitution.set_multiplicity(multiplicity)
                defects.append(substitution)
    return defects   


def create_vacancies(structure,elements=None,supercell_size=None, **kwargs):
    """
    Create Vacancy objects for each non-equivalent site,
    staring from a bulk structure (unit cell or supercell).

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
    defects : (list)
        List of Vacancy objects
    """
    defects = []
    bulk_structure = structure.copy()
    if supercell_size:
        bulk_structure.make_supercell(supercell_size)
    if not elements:
        elements = [el.symbol for el in bulk_structure.composition.elements]
        
    sym_struct = SpacegroupAnalyzer(structure=bulk_structure,**kwargs).get_symmetrized_structure()
    for el in bulk_structure.composition.elements:
        for site_group in sym_struct.equivalent_sites:
            site = site_group[0]
            if el.symbol in elements:
                if site.specie == el:
                    vacancy = Vacancy(
                                    defect_site=site,
                                    bulk_structure=bulk_structure)
                    multiplicity = len(sym_struct.find_equivalent_sites(site))
                    vacancy.set_multiplicity(multiplicity)
                    defects.append(vacancy)
    return defects     