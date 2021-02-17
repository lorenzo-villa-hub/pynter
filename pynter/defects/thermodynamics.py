#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:48:55 2021

@author: villa
"""

from pynter.defects.analysis import DefectsAnalysis
from pynter.phase_diagram.experimental import ChempotExperimental



class PartialPressureAnalysis:
    """
    Class that handles the analysis of the oxygen partial pressure dependency.
    """
    
    def __init__(self,defects_analysis,phase_diagram,target_comp,bulk_dos,temperature,frozen_defect_concentrations=None,external_defects=[]):
        """
        Parameters
        ----------
        defects_analysis :
            DefectsAnalysis object.
        phase_diagram : 
            Pymatgen PhaseDiagram object.
        target_comp : 
            Pymatgen Composition object.
        bulk_dos : 
            Pymatgen Dos object.
        temperature : (float), optional
            Temperature in Kelvin. The default is 300K.
        """
        self.da = defects_analysis
        self.pd = phase_diagram
        self.target_comp = target_comp
        self.bulk_dos = bulk_dos
        self.temperature = temperature
        self.frozen_defect_concentrations = frozen_defect_concentrations
        self.external_defects = external_defects
    
    
    def get_concentrations(self,pressure_range=(-20,10),concentrations_output='all',npoints=50):
        """
        Calculate defect and carrier concentrations at different oxygen partial pressure values

        Parameters
        ----------
        pressure_range : (tuple), optional
            Exponential range in which to evaluate the partial pressure. The default is from 1e-20 to 1e10.
        concentrations_output : (str), optional
            Type of output for defect concentrations:
                "all": The output is the concentration of every defect entry.
                "stable": The output is the concentration of the stable charge for every defect at each fermi level point.
                "total": The output is the sum of the concentration in every charge for each specie.
                The default is 'all'.

        Returns
        -------
        partial_pressures : (list)
            List of partial pressure values.
        defect_concentrations : (list)
            If the output is set to "all" is a list of list of dictionaries with "name", "charge", "conc" as keys. 
            If the output is "all" is a list of dictionaries with names as keys and conc as values. 
        carrier_concentrations : (list)
            List of tuples with intrinsic carriers concentrations (holes,electrons).
        """
        
        res = ChempotExperimental().chempots_partial_pressure_range(self.pd,self.target_comp,
                                                                  self.temperature,pressure_range=pressure_range,npoints=npoints)
        partial_pressures = list(res.keys())
        defect_concentrations = []
        carrier_concentrations = []
        dos = self.bulk_dos
        T = self.temperature
        frozen_df = self.frozen_defect_concentrations
        ext_df = self.external_defects
        for r,mu in res.items():
            if frozen_df or ext_df:
                mue = self.da.non_equilibrium_fermi_level(frozen_df,mu,dos,ext_df,temperature=T)
            else:
                mue = self.da.equilibrium_fermi_level(mu,dos,temperature=T)
            if concentrations_output == 'all':
                conc = self.da.defect_concentrations(mu,temperature=T,fermi_level=mue)
            elif concentrations_output == 'total':
                conc = self.da.defect_concentrations_total(mu,temperature=T,fermi_level=mue)
            elif concentrations_output == 'stable':
                conc = self.da.defect_concentrations_stable_charges(mu,temperature=T,fermi_level=mue)
            else:
                raise ValueError('concentrations_output must be chosen between "all", "total", "stable"') 
            carriers = self.da.carrier_concentrations(dos,temperature=T,fermi_level=mue)
            defect_concentrations.append(conc)
            carrier_concentrations.append(carriers)
            
        return partial_pressures, defect_concentrations, carrier_concentrations
    
    
    def get_fermi_levels(self,pressure_range=(-20,10),npoints=50):
        """
        Calculate defect and carrier concentrations at different oxygen partial pressure values

        Parameters
        ----------
        pressure_range : (tuple), optional
            Exponential range in which to evaluate the partial pressure. The default is from 1e-20 to 1e10.
            
        Returns
        -------
        partial_pressures : (list)
            List of partial pressure values.
        fermi_levels : (list)
            List of Fermi level values
        """
        res = ChempotExperimental().chempots_partial_pressure_range(self.pd,self.target_comp,
                                                                  self.temperature,pressure_range=pressure_range,npoints=npoints)
        partial_pressures = list(res.keys())
        fermi_levels = []
        dos = self.bulk_dos
        T = self.temperature
        frozen_df = self.frozen_defect_concentrations
        ext_df = self.external_defects
        for r,mu in res.items():
            if frozen_df or ext_df:
                mue = self.da.non_equilibrium_fermi_level(frozen_df,mu,dos,ext_df,temperature=T)
            else:
                mue = self.da.equilibrium_fermi_level(mu,dos,temperature=T)
            fermi_levels.append(mue)
            
        return partial_pressures, fermi_levels
    

    





