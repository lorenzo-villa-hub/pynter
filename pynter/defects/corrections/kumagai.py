
import warnings
import os.path as op
import importlib
import matplotlib.pyplot as plt
import math

from pynter.tools.structure import is_site_in_structure_coords
from pynter.defects.structure import defect_finder

from pymatgen.io.vasp.outputs import Outcar, Vasprun

from pymatgen.analysis.defects.corrections.kumagai import _check_import_pydefect 
from pymatgen.analysis.defects.utils import CorrectionResult, get_zfile

from pathlib import Path



    
def get_kumagai_correction(
                        defect_path,
                        bulk_path,
                        charge,
                        dielectric_tensor,
                        initial_structure=False,
                        get_correction_data=True,
                        get_plot=True,
                        **kwargs):
    """
    Compute Kumagai corrections (extended FNV scheme) from VASP calculation paths.

    Parameters
    ----------
    defect_path : (str)
        Path of defect calculation with vasprun.xml and OUTCAR files.
    bulk_path : (str)
        Path of bulk calculation with vasprun.xml and OUTCAR files.
    charge : (int)
        Charge of defect calculation
    dielectric_tensor : (int,float 3x1 array or 3x3 array)
        Dielectric tensor (or constant). Types accepted are int,float 3x1 array or 3x3 array.
    initial_structure : (bool)
        Use initial structure of defect calculation for correction computation.
        Useful to compute correction on unrelaxed structure.
    get_correction_data : (bool)
        Return pymatgen's CorrectionResult object. If False only the correction value is returned.
    get_plot : (bool), optional
        Get Matplotlib object with plot. The default is False. 
    kwargs : (dict)
        Kwargs to pass to pydefect `make_efnv_correction`
    
    Returns
    -------
    correction, ax
        - CorrectionResult object if get_correction_data is set to True, else just the float with correction value,
        - matplotlib axis object  
    """
    defect_structure = _get_structure_with_pot_pmg(defect_path,initial_structure=initial_structure)
    bulk_structure = _get_structure_with_pot_pmg(bulk_path,initial_structure=False)
    correction = get_kumagai_correction_from_structures(
                                            defect_structure_with_potentials=defect_structure,
                                            bulk_structure_with_potentials=bulk_structure,
                                            charge=charge,
                                            dielectric_tensor=dielectric_tensor,
                                            **kwargs)
    
    corr = correction if get_correction_data else correction.correction_energy
    if get_plot:
        from pydefect.corrections.site_potential_plotter import SitePotentialMplPlotter
        SitePotentialMplPlotter.from_efnv_corr(
                                            title=defect_structure.composition,
                                            efnv_correction=correction.metadata['efnv_corr']
                                            ).construct_plot()
        ax = plt.gca()
        return corr, ax
    else:
        return corr 


def get_kumagai_correction_from_jobs(
                                    job_defect,
                                    job_bulk,
                                    dielectric_tensor,
                                    initial_structure=False,
                                    get_correction_data=True,
                                    get_plot=True,
                                    **kwargs):
    """
    Get Kumagai corrections from VaspJob objects.

    Parameters
    ----------
    job_defect : (VaspJob)
        Defect calculation.
    job_bulk : (VaspJob)
        Bulk calculation.
    dielectric_tensor : (int,float 3x1 array or 3x3 array)
        Dielectric tensor (or constant). Types accepted are int,float 3x1 array or 3x3 array.
    initial_structure : (bool)
        Use initial structure of defect calculation for correction computation.
        Useful to compute correction on unrelaxed structure.
    get_correction_data : (bool)
        Return pymatgen's CorrectionResult object. If False only the correction value is returned.
    get_plot : (bool), optional
        Get Matplotlib object with plot. The default is False. 
    kwargs : (dict)
        Kwargs to pass to pydefect `make_efnv_correction`
    
    Returns
    -------
    correction, ax
        - CorrectionResult object if get_correction_data is set to True, else just the float with correction value,
        - matplotlib axis object  
    """
    defect_path, bulk_path = job_defect.path, job_bulk.path
    charge = job_defect.charge

    corr = get_kumagai_correction(
                        defect_path,
                        bulk_path,
                        charge,
                        dielectric_tensor,
                        initial_structure=False,
                        get_plot=False,
                        **kwargs)    
    return corr


def get_kumagai_correction_from_structures(
                                        defect_structure_with_potentials,
                                        bulk_structure_with_potentials,
                                        charge,
                                        dielectric_tensor,
                                        **kwargs):
    """
    Get Kumagai (extended FNV) correction from Structure objects with 
    site potentials stored into site_properties. Not recommended to use
    directly, refer to get_kumagai_correction.
    """
    dielectric_tensor = _convert_dielectric_tensor(dielectric_tensor)
    correction =  _get_efnv_correction_pmg_fixed(
                            charge=charge,
                            defect_structure=defect_structure_with_potentials,
                            bulk_structure=bulk_structure_with_potentials,
                            dielectric_tensor=dielectric_tensor,
                            **kwargs)
    return correction




def _get_efnv_correction_pmg_fixed(
    charge,
    defect_structure,
    bulk_structure,
    dielectric_tensor,
    **kwargs):
    """
    Fixed code from pymatgen.analysis.defects.corrections.kumagai.get_efnv_correction
    defect_potentials and site_potentials variables were inverted.
    
    Returns the Kumagai/Oba EFNV correction for a given defect.

    Args:
        charge: Charge of the defect.
        defect_structure: Defect structure.
        bulk_structure: Bulk structure.
        dielectric_tensor: Dielectric tensor.
        **kwargs: Keyword arguments to pass to `make_efnv_correction`.
    """
    from pydefect.analyzer.calc_results import CalcResults
    from pydefect.cli.vasp.make_efnv_correction import make_efnv_correction

    # ensure that the structures have the "potential" site property
    defect_potentials = [site.properties["potential"] for site in defect_structure]
    bulk_potentials = [site.properties["potential"] for site in bulk_structure]

    defect_calc_results = CalcResults(
        structure=defect_structure,
        energy=math.inf,
        magnetization=math.inf,
        potentials=defect_potentials,
    )
    bulk_calc_results = CalcResults(
        structure=bulk_structure,
        energy=math.inf,
        magnetization=math.inf,
        potentials=bulk_potentials,
    )

    efnv_corr = make_efnv_correction(
        charge=charge,
        calc_results=defect_calc_results,
        perfect_calc_results=bulk_calc_results,
        dielectric_tensor=dielectric_tensor,
        **kwargs,
    )

    return CorrectionResult(
        correction_energy=efnv_corr.correction_energy,
        metadata={"efnv_corr": efnv_corr},
    )


def _get_structure_with_pot_pmg(directory,initial_structure=False):
    """
    Modified function from pymatgen.analysis.defects-corrections.kumagai 
    to allow for initial structure import

    Reads vasprun.xml and OUTCAR files in a directory.

    Args:
        directory: Directory containing vasprun.xml and OUTCAR files.

    Returns:
        Structure with "potential" site property.
    """
    _check_import_pydefect()
    from pydefect.analyzer.calc_results import CalcResults

    d_ = Path(directory)
    f_vasprun = get_zfile(d_, "vasprun.xml")
    f_outcar = get_zfile(d_, "OUTCAR")
    vasprun = Vasprun(f_vasprun,parse_dos=False,parse_potcar_file=False,parse_eigen=False)
    outcar = Outcar(f_outcar)

    structure = vasprun.structures[0] if initial_structure else vasprun.final_structure

    calc = CalcResults(
        structure=structure,
        energy=outcar.final_energy,
        magnetization=outcar.total_mag or 0.0,
        potentials=[-p for p in outcar.electrostatic_potential],
        electronic_conv=vasprun.converged_electronic,
        ionic_conv=vasprun.converged_ionic,
    )

    return calc.structure.copy(site_properties={"potential": calc.potentials})


def _convert_dielectric_tensor(dielectric):
    import numpy as np
    if not type(dielectric) in [float,int]:
        dielectric = np.array(dielectric)
        if dielectric.shape == (3,):
            dielectric = np.diag(dielectric)
        elif dielectric.shape != (3, 3):
            raise ValueError("Dielectric tensor can be int/float, a 3x1 array with diagonal components of dielectric tensor, or 3x3 matrix")
    else:
        dielectric = np.eye(3) * dielectric

    return dielectric



# ########################## CORRECTIONS WITH OLD pymatgen defects #########################################

# from ._pmg.pmg_defects_core import Vacancy , DefectEntry, Interstitial, Substitution
# from ._pmg.pmg_defects_corrections import KumagaiCorrection


# def get_kumagai_correction_old_pmg(
#                         defect_path,
#                         bulk_path,
#                         charge,
#                         dielectric_tensor,
#                         sampling_radius=None,
#                         gamma=None,
#                         get_plot=False,
#                         initial_structure=True,
#                         **kwargs):
#     vasprun_defect = Vasprun(op.join(defect_path,'vasprun.xml'),parse_potcar_file=False,parse_dos=False,parse_eigen=False)
#     vasprun_bulk = Vasprun(op.join(defect_path,'vasprun.xml'),parse_potcar_file=False,parse_dos=False,parse_eigen=False)
#     if initial_structure:
#         structure_defect = vasprun_defect.structures[0]
#     else:
#         structure_defect = vasprun_defect.final_structure
#     structure_bulk = vasprun_bulk.final_structure
#     path_to_defect_outcar = op.join(defect_path,'OUTCAR')
#     path_to_bulk_outcar = op.join(bulk_path,'OUTCAR')

#     corr = get_kumagai_correction_from_old_pmg(
#                                             structure_defect=structure_defect,
#                                             structure_bulk=structure_bulk,
#                                             path_to_defect_outcar=path_to_defect_outcar,
#                                             path_to_bulk_outcar=path_to_bulk_outcar,
#                                             dielectric_tensor=dielectric_tensor,
#                                             charge=charge,
#                                             sampling_radius=sampling_radius,
#                                             gamma=gamma,
#                                             get_plot=get_plot,
#                                             **kwargs)
    
#     return corr
    

# def get_kumagai_correction_from_old_pmg(
#                         structure_defect,
#                         structure_bulk,
#                         path_to_defect_outcar,
#                         path_to_bulk_outcar,
#                         dielectric_tensor,
#                         charge,
#                         defect_type=None,
#                         defect_specie=None,
#                         defect_site=None,
#                         sampling_radius=None,
#                         gamma=None,
#                         get_plot=False,
#                         **kwargs):
#     """
#     Get Kumagai correction with Pymatgen.

#     Parameters
#     ----------
#     structure_defect : (Structure)
#         Structure of defect.
#     structure_bulk : (Structure)
#         Bulk structure.
#     path_to_defect_outcar : (str)
#         Path to OUTCAR of defect calculation.
#     path_to_bulk_outcar : (str)
#         Path to OUTCAR of pure calculation.
#     dielectric_tensor : (array or float)
#         Dielectric tensor, if is a float a diagonal matrix is constructed.
#     charge : (int or float)
#         Charge of the defect.
#     defect_type : (str), optional
#         Type of defect ('Vacancy','Interstitial' or 'Substitution')
#         If None it's determined with defect_finder. The default is None.
#     defect_specie : (str), optional
#         Symbol of the defect specie.
#         If None it's determined with defect_finder. The default is None.
#     defect_site : (Site), optional
#         Site of defect. If None it's determined with defect_finder. The default is None.
#     sampling_radius (float): radius (in Angstrom) which sites must be outside
#         of to be included in the correction. Publication by Kumagai advises to
#         use Wigner-Seitz radius of defect supercell, so this is default value.
#     gamma (float): convergence parameter for gamma function.
#                     Code will automatically determine this if set to None.
#     tol : (float)
#         Tolerance for comparing sites and defect_finder function. The default is 1e-03.
#     get_plot : (bool), optional
#         Get Matplotlib object with plot. The default is False.

#     Returns
#     -------
#     corr : (dict or tuple)
#         Dictionary with corrections, if get_plot is True a tuple with dict and plt object is returned.
#     """
#     if 'max_number_of_defects' not in kwargs.keys():
#         kwargs['max_number_of_defects'] = 1
#     elif kwargs['max_number_of_defects'] != 1:
#         kwargs['max_number_of_defects'] = 1
#         warnings.warn('Corrections are only defined for single defects')
    
#     if not defect_site and not defect_type and not defect_specie:
#         df_found = defect_finder(structure_defect, structure_bulk, **kwargs)
#         defect_site, defect_type = df_found.site, df_found.type
#         defect_specie = defect_site.specie.symbol
   
#     is_site_kwargs = {}
#     if 'tol' in kwargs:
#         is_site_kwargs['tol'] = kwargs['tol']
#     site_matching_indices = []
#     for site in structure_defect:
#         site_in_str ,index_bulk = is_site_in_structure_coords(site, structure_bulk, **is_site_kwargs)
#         if site_in_str:
#             site_matching_indices.append([index_bulk,structure_defect.index(site)])
#         else:
#             print(f'Warning in Kumagai corrections: Site {site} is not in bulk structure')
    
#     bulk_atomic_site_averages = Outcar(path_to_bulk_outcar).read_avg_core_poten()[-1]
#     defect_atomic_site_averages = Outcar(path_to_defect_outcar).read_avg_core_poten()[0]
#     defect_frac_sc_coords = defect_site.frac_coords
#     initial_defect_structure = structure_defect
    
#     parameters = {}
#     parameters['bulk_atomic_site_averages'] = bulk_atomic_site_averages
#     parameters['defect_atomic_site_averages'] = defect_atomic_site_averages
#     parameters['site_matching_indices'] = site_matching_indices
#     parameters['initial_defect_structure'] = initial_defect_structure
#     parameters['defect_frac_sc_coords'] = defect_frac_sc_coords
    
#     module = importlib.import_module("pynter.defects.corrections._pmg.pmg_defects_core")
#     defect_class = getattr(module,defect_type)
#     defect = defect_class(structure_bulk, defect_site, charge=charge, multiplicity=1)
#     defect_entry = DefectEntry(defect,None,corrections=None,parameters=parameters)

#     kumagai = KumagaiCorrection(dielectric_tensor,sampling_radius,gamma)
#     kumagai_corrections = kumagai.get_correction(defect_entry)
    
#     if get_plot:
#         plt = kumagai.plot()
#         return kumagai_corrections , plt
#     else:    
#         return kumagai_corrections