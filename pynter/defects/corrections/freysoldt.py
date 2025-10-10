

import os.path as op
import importlib

from pymatgen.io.vasp.outputs import Locpot
from pymatgen.analysis.defects.corrections.freysoldt import get_freysoldt_correction, plot_plnr_avg

from ..structure import defect_finder


def get_freysoldt_correction_from_locpot(
                                    charge,
                                    dielectric_constant,
                                    defect_path_locpot,
                                    bulk_path_locpot,
                                    get_correction_data=True,
                                    get_plot=True,
                                    plot_axis_index=0,
                                    finder_kwargs={},
                                    **kwargs):
    """
    Get Freysoldt correction from defect and bulk LOCPOT files, using "pymatgen-analysis-defects".

    Parameters
    ----------
    charge : (int)
        Charge of the defect
    dielectric_constant : (float)
        Dielectric constant of the material
    defect_path_locpot : (str)
        Defect LOCPOT file.
    bulk_path_locpot : (str)
        Bulk LOCPOT file.
    get_correction_data : (bool)
        Return pymatgen's CorrectionResult object. If False only the correction value is returned.
    get_plot : (bool), optional
        Get Matplotlib object with plot. The default is False.
    plot_axis_index : (int), optional
        Axis for planar average plot. The default is 0.
    finder_kwargs : (dict)
        Kwargs to pass to defect_finder.
    kwargs : (dict)
        Kwargs to pass to pymatgen's get_freysoldt_correction.

    Returns
    -------
    correction, ax
        - CorrectionResult object if get_correction_data is set to True, else just the float with correction value,
        - matplotlib axis object
    """
    q = charge
    defect_locpot = Locpot.from_file(defect_path_locpot)
    bulk_locpot = Locpot.from_file(bulk_path_locpot)
    if finder_kwargs:
        if 'verbose' not in finder_kwargs.keys():
            finder_kwargs['verbose'] = True
    defect = defect_finder(defect_locpot.structure,bulk_locpot.structure,**finder_kwargs)
    defect_frac_coords = defect.site.frac_coords
    lattice = defect_locpot.structure.lattice
    correction = get_freysoldt_correction(q=q,dielectric=dielectric_constant,defect_locpot=defect_locpot,
                                    bulk_locpot=bulk_locpot,defect_frac_coords=defect_frac_coords,
                                    lattice=lattice,**kwargs)
    
    corr = correction if get_correction_data else correction.correction_energy
    if get_plot:
        ax = plot_plnr_avg(correction.metadata['plot_data'][plot_axis_index])
        return corr, ax
    else:
        return corr 
    


def get_freysoldt_correction_from_jobs(
                                    job_defect,
                                    job_bulk,
                                    dielectric_constant,
                                    get_correction_data=True,
                                    get_plot=False,
                                    plot_axis_index=0,
                                    finder_kwargs={},
                                    **kwargs):
    """
    Get Freysoldt correction from VaspJob objects using the new "pymatgen-analysis-defects".

    Parameters
    ----------
    job_defect : (VaspJob)
        Defect calculation.
    job_bulk : (VaspJob)
        Bulk calculation.
    dielectric_constant : (float)
        Dielectric constant.
    get_correction_data : (bool)
        Return pymatgen's CorrectionResult object. If False only the correction value is returned.
    get_plot : (bool), optional
        Get Matplotlib object with plot. The default is False.
    plot_axis_index : (int), optional
        Axis for planar average plot. The default is 0.
    finder_kwargs : (dict)
        Kwargs to pass to defect_finder.
    kwargs : (dict)
        Kwargs to pass to pymatgen's get_freysoldt_correction.

    Returns
    -------
    correction, ax
        - CorrectionResult object if get_correction_data is set to True, else just the float with correction value,
        - matplotlib axis object
    """
    q = job_defect.charge
    defect_locpot = Locpot.from_file(op.join(job_defect.path,'LOCPOT'))
    bulk_locpot = Locpot.from_file(op.join(job_bulk.path,'LOCPOT'))
    return get_freysoldt_correction_from_locpot(
                                                charge=q,
                                                dielectric_constant=dielectric_constant,
                                                defect_path_locpot=defect_locpot,
                                                bulk_path_locpot=bulk_locpot,
                                                get_correction_data=get_correction_data,
                                                get_plot=get_plot,
                                                plot_axis_index=plot_axis_index,
                                                finder_kwargs=finder_kwargs,
                                                **kwargs)