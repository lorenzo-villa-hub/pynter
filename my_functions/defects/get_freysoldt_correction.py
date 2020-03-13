#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:41:35 2019

@author: villa
"""
##################################################################################################

def get_freysoldt_correction(defect_type, defect_specie, path_to_defect_locpot,path_to_pure_locpot,charge,
                             dielectric_constant,defect_site_coordinates,energy_cutoff=500):
    
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
            
        Returns:
            Freysoldt corrections values as a dictionary 
            '''
    
    import numpy as np
    from pymatgen.core.sites import PeriodicSite
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.io.vasp.outputs import Vasprun, Locpot, VolumetricData
    from pymatgen.analysis.defects.core import Vacancy , DefectEntry, Interstitial
    from pymatgen.analysis.defects.corrections import FreysoldtCorrection
    from pymatgen.analysis.defects.generators import SubstitutionGenerator

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
    
#    # finding defect site
#    for site in locpot_pure.structure.sites:
#        if np.array_equiv(np.around(site.frac_coords,3),np.around(defect_site_coordinates,3)):
#            defect_site = site
    
    #setting defect site
    defect_site = PeriodicSite(defect_specie, coords=defect_site_coordinates, lattice = locpot_pure.structure.lattice)
    
    if defect_type == 'vacancy':
        defect = Vacancy(locpot_pure.structure, defect_site, charge=charge, multiplicity=None)
    if defect_type == 'interstitial':
        defect = Interstitial(locpot_pure.structure, defect_site, charge=charge, multiplicity=None)
    
    defect_entry = DefectEntry(defect,None,corrections=None,parameters=parameters)
    
    freysoldt_class = FreysoldtCorrection(dielectric_constant,energy_cutoff=energy_cutoff)
    
    freysoldt_corrections = freysoldt_class.get_correction(defect_entry)
  
  # plot  
  #  freysoldt_class.plot(1)
    
    return freysoldt_corrections

##############################################################################
    
