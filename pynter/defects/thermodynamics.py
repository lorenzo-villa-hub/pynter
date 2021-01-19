#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:48:55 2021

@author: villa
"""

from pynter.defects.analysis import DefectsAnalysis
from pynter.phase_diagram.experimental import ChempotExperimental
import matplotlib
import matplotlib.pyplot as plt
from pymatgen.util.plotting import pretty_plot


class PartialPressureAnalysis:
    """
    Class that handles the analysis of the oxygen partial pressure dependency.
    """
    
    def __init__(self,defects_analysis,phase_diagram,target_comp,bulk_dos,temperature=300,frozen_defect_concentrations=None,external_defects=[]):
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
    
    
    def get_concentrations(self,pressure_range=(-20,10),concentrations_output='all'):
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
        
        res = ChempotExperimental().oxygen_partial_pressure_range(self.pd,self.target_comp,
                                                                  self.temperature,p_range=pressure_range)
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
    
    
    def get_fermi_levels(self,pressure_range=(-20,10)):
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
        res = ChempotExperimental().oxygen_partial_pressure_range(self.pd,self.target_comp,
                                                                  self.temperature,p_range=pressure_range)
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
    
    
    def plot_concentrations(self,partial_pressures,defect_concentrations,carrier_concentrations,
                            defect_indexes=None,concentrations_output='all',size=(12,8),xlim=(1e-20,1e08),ylim=None):
        """
        Plot defect and carrier concentrations in a range of oxygen partial pressure.

        Parameters
        ----------
        defect_indexes : (list), optional
            List of indexes of the entry that need to be included in the plot. If None 
            all defect entries are included. The default is None.
        pressure_range : (tuple), optional
            Exponential range in which to evaluate the partial pressure. The default is from 1e-20 to 1e10.
        concentrations_output : (str), optional
            Type of output for defect concentrations:
                "all": The output is the concentration of every defect entry.
                "stable": The output is the concentration of the stable charge for every defect at each fermi level point.
                "total": The output is the sum of the concentration in every charge for each specie.
                The default is 'all'.
        size : (tuple), optional
            Size of the matplotlib figure. The default is (12,8).
        xlim : (tuple), optional
            Range of x-axis. The default is (1e-20,1e08).
        ylim : (tuple), optional
            Range of y-axis. The default is None.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        p,dc,cc = partial_pressures,defect_concentrations,carrier_concentrations
        matplotlib.rcParams.update({'font.size': 22})
        if concentrations_output == 'all' or concentrations_output == 'stable':
            plt = self._plot_conc(p,dc,cc,defect_indexes,concentrations_output,size)
        elif concentrations_output == 'total':
            plt = self._plot_conc_total(p,dc,cc,size)
            
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        plt.xlabel('Oxygen partial pressure (atm)')
        plt.ylabel('Concentrations(cm$^{-3})$')
        plt.legend()
        plt.grid()
        
        return plt
        
   
    def plot_fermi_level(self,partial_pressures, fermi_levels, new_figure=True,label=None,size=(12,8),xlim=(1e-20,1e08),ylim=None):
        """
        Plot Fermi level as a function of the oxygen partial pressure.

        Parameters
        ----------
        new_figure : (bool), optional
            Initialize a new matplotlib figure. The default is True.
        label : (str), optional
            Label for the data. The default is None.
        pressure_range : (tuple), optional
            Exponential range in which to evaluate the partial pressure. The default is from 1e-20 to 1e10.
        size : (tuple), optional
            Size of the matplotlib figure. The default is (12,8).
        xlim : (tuple), optional
            Range of x-axis. The default is (1e-20,1e08).
        ylim : (tuple), optional
            Range of y-axis. The default is None.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        if not label:
            label = '$\mu_{e}$'
        matplotlib.rcParams.update({'font.size': 22})
        if new_figure:
            plt.figure(figsize=(size))
        p,mue = partial_pressures, fermi_levels
        plt.plot(p,mue,linewidth=4,label=label)
        plt.xscale('log')
        plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        plt.xlabel('Oxygen partial pressure (atm)')
        plt.ylabel('Electron chemical potential (eV)')
        plt.legend()
        if new_figure:
            plt.grid()
        
        return plt
    
    
    def _plot_conc(self,partial_pressures,defect_concentrations,carrier_concentrations,defect_indexes,
                   concentrations_output,size):
        
        plt.figure(figsize=size)
        filter_defects = True if defect_indexes else False
        p = partial_pressures
        dc = defect_concentrations
        h = [cr[0] for cr in carrier_concentrations] 
        n = [cr[1] for cr in carrier_concentrations] 
        for i in range(0,len(defect_concentrations[0])):
            if filter_defects is False or (defect_indexes is not None and i in defect_indexes):
                conc = [c[i]['conc'] for c in defect_concentrations]
                label_txt = self.da._get_formatted_legend(dc[0][i]['name'])
                label_txt = self._format_legend_with_charge(label_txt,dc[0][i]['charge'])
                plt.plot(p,conc,label=label_txt,linewidth=4)
        plt.plot(p,h,label='$n_{h}$',linestyle='--',color='r',linewidth=4)
        plt.plot(p,n,label='$n_{e}$',linestyle='--',color='b',linewidth=4)
        return plt
    
    
    def _plot_conc_total(self,partial_pressures,defect_concentrations,carrier_concentrations,size):
        
        plt.figure(figsize=size)
        p = partial_pressures
        h = [cr[0] for cr in carrier_concentrations] 
        n = [cr[1] for cr in carrier_concentrations] 
        for name in defect_concentrations[0]:
            conc = [c[name] for c in defect_concentrations]
            label_txt = self.da._get_formatted_legend(name)
            plt.plot(p,conc,label=label_txt,linewidth=4)
        plt.plot(p,h,label='$n_{h}$',linestyle='--',color='r',linewidth=4)
        plt.plot(p,n,label='$n_{e}$',linestyle='--',color='b',linewidth=4)
        return plt
    
    
    def _format_legend_with_charge(self,label,charge):
        
        mod_label = label[:-1]
        if charge < 0:
            for i in range(0,abs(charge)):
                if i == 0:
                    mod_label = mod_label + "^{"
                mod_label = mod_label + "°"
            mod_label = mod_label + "}"
        elif charge == 0:
            mod_label = mod_label + "^{x}"
        elif charge > 0:
            for i in range(0,charge):
                if i == 0:
                    mod_label = mod_label + "^{"
                mod_label = mod_label + "'"
            mod_label = mod_label + "}"
        
        mod_label = mod_label + "$"
        return mod_label
    

    
    





