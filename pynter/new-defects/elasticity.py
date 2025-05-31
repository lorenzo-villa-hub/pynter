#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 14:35:23 2021

@author: villa
"""

import numpy as np
import pandas as pd
import matplotlib

class Stresses:
    
    def __init__(self,stress_bulk):
        """
        Class to handle elastic analysis from defect calculations. Necessary data must be present for 
        defect entries in dictionary form in "data" attribute.
        Valid keys are:
            "stress" : To compute elastic dipole tensor and relaxation volume starting from stress tensor 
                    of defect calculation.
            "elastic_corrections" : dict with corrections for charged defects ({"corr_name":corr_value}).
                    This values will be added to the elastic dipole tensor.

        Parameters
        ----------
        stress_bulk : (np.array)
            Stress tensor of bulk calculation.
        """
        self.stress_bulk = stress_bulk
        
        
    def get_elastic_dipole_tensor(self,entry,add_corrections=True,add_to_data=False):
        """
        Calculate elastic dipole tensor as P_ij = -Omega * sigma_ij, where Omega is the 
        bulk volume and sigma_ij is the residual stress tensor.

        Parameters
        ----------
        entry : (DefectEntry)
            DefectEntry object. Stress needs to be present in the data
            dictionary as a numpy array with key "stress".
        add_corrections : (bool)
            Add correction terms from "elastic_corrections" dict (if key is present in dict). The default is True.
 
        Returns
        -------
        dipole_tensor : (np.array)
            Elastic dipole tensor in eV.
        """
        bulk_volume = entry.bulk_structure.lattice.volume
        res_stress = self.get_residual_stress_tensor(entry)
        dipole_tensor = -1*bulk_volume*res_stress #sign is inverted with respect to VASP output
        if 'elastic_corrections' in entry.data.keys() and add_corrections:
            dipole_tensor += sum([v for k,v in entry.data['elastic_corrections'].items()]) 
        if add_to_data:
            entry.data['elastic_dipole_tensor'] = dipole_tensor
        return dipole_tensor


    def get_relaxation_volume(self,entry,bulk_modulus,add_corrections=True,add_to_data=True):
        """
        Calculate relaxation volume from stresses. Stresses data needs to be in numpy.array format and present 
        in the "data" dictionary with relative "stress" key.

        Parameters
        ----------
        entry : (DefectEntry)
            DefectEntry object. Stress needs to be present in the data dictionary 
            as a numpy array with key "stress".  
        bulk_modulus : (float)
            Bulk modulus in GPa.
        add_corrections : (bool)
            Add correction terms from "elastic_corrections" dict (if key is present in dict). The default is True.
        add_to_data : (bool)
            Add relaxation volume value to entry.data dict with key "relaxation_volume". The default is True.
            
        Returns
        -------
        rel_volume : (float)
            Relaxation volume in A°^3.
        """
        bulk_modulus = bulk_modulus*10 # from GPa to kbar
        dipole_tensor = self.get_elastic_dipole_tensor(entry)
        pressure = np.trace(dipole_tensor)/3
        rel_volume = pressure/bulk_modulus
        if add_to_data:
            entry.data['relaxation_volume'] = rel_volume
        return rel_volume

        
    def get_residual_stress_tensor(self,entry): 
        """
        Calculate residual stress tensor by subtracting bulk stress tensor.

        Parameters
        ----------
        entry : (DefectEntry)
            DefectEntry object. Stress needs to be present in the data dictionary 
            as a numpy array with key "stress".

        Returns
        -------
        stress : (np.array)
            Residual stress matrix.
        """ 
        if 'stress' in entry.data.keys():
            stress_def = entry.data['stress']
            stress = np.array(stress_def) - np.array(self.stress_bulk) # residual stress -> defect - bulk
            return stress
        else:
            raise UserWarning('Stress needs to be present in the data dictionary as a numpy array with key "stress".')
            return
    
        
    def plot_relaxation_volumes(self,entries,bulk_modulus,check_in_data=True,add_corrections=True,
                                add_to_data=True,get_dataframe=False,**kwargs):
        """
        Plot relaxation volumes for list of defect entries with pandas Series.

        Parameters
        ----------
        entries : (list)
            List of DefectEntry objects.
        bulk_modulus : (float)
            Bulk modulus in GPa.
        check_in_data : (bool), optional
            Check weather relaxation volume is already present in entry.data['relaxation_volume']. The default is True.
        add_corrections : (bool)
            Add correction terms from "elastic_corrections" dict (if key is present in dict). The default is True.
        add_to_data : (bool)
            Add relaxation volume value to entry.data dict with key "relaxation_volume". The default is True.
        get_dataframe : (bool), optional
            Return tuple with pandas data: (plt,dataframe). If False only the matplotlib object is returned.
            The default is False.
        **kwargs : (dict)
            Kwargs for df.plot method.

        Returns
        -------
        (matplotlib object)
            If get_dataframe is True return Tuple with matplotlib object and pandas data: (plt,dataframe). 
            If False only the matplotlib object is returned.
        """
        rel_volumes = {}
        for e in entries:
            if check_in_data and 'relaxation_volume' in e.data.keys():
                rel_volumes[e.symbol] = e.data['relaxation_volume']                  
            else:
                rel_volumes[e.symbol] = self.get_relaxation_volume(e,bulk_modulus,add_corrections,add_to_data)
                
        df = pd.Series(rel_volumes)
        if not kwargs:
            kwargs = {'grid':True,'figsize':(12,8)} #default values
        fontsize = 15*kwargs['figsize'][0]/12 
        matplotlib.rcParams.update({'font.size':fontsize})
        plt = df.plot.bar(rot=0,ylabel='Relaxation Volume ($\AA^{3}$)',**kwargs).get_figure()
        
        if get_dataframe:
            return plt,df
        else:
            return plt
        
    
    
    
    
    