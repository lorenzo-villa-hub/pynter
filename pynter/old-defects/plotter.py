#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 15:24:29 2021

@author: villa
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from pynter.defects.defects import format_legend_with_charge_number, get_defect_name_from_string


class ConcBarPlotter:

    def __init__(self,concentrations,format_names=True):
        """
        Class to visualize concentrations dictionaries with pd.DataFrame and 
        pd.Series (fortotal concentrations).

        Parameters
        ----------
        concentrations : (list or dict)
            Concentrations can be in different formats: list of dict for "all" and 
            "stable_charges" and dict for "total".
        format_names : (bool), optional
            Format names with latex symbols. The default is True.
        """
        # to be fixed
        conc_dict = []
        for c in concentrations:
            d = {'charge':c.charge,'conc':c.conc,'stable':c.stable}
            if format_names:
                name = c.name.symbol
            else:
                name = c.name.fullname
            d['name'] = format_legend_with_charge_number(name,c.charge)
            conc_dict.append(d)
        
        conc_total_dict = {}
        for dn,conc in concentrations.total.items():
            if format_names:
                name = dn.symbol
            else:
                name = dn.fullname
            conc_total_dict[name] = conc

        self.conc = conc_dict
        self.conc_total = conc_total_dict
        self.format = format_names
        self.df = pd.DataFrame(self.conc)
        self.series = pd.Series(self.conc_total,name='Total Concentrations')
    

    def __print__(self):
        return self.df.__print__()
    
    def __repr__(self):
        return self.df.__repr__()

    def copy(self):
        return self.df.copy()
    
    
    def plot_bar(self,conc_range=(1e13,1e40),ylim=None,total=True,ylabel_fontsize=15,**kwargs):
        """
        Bar plot of concentrations with pd.DataFrame.plot

        Parameters
        ----------
        conc_range : (tuple), optional
            Range of concentrations to include in df. The default is (1e13,1e40).
        ylim : (tuple), optional
            Limit of y-axis in plot. If None conc_range is used. The default is None.
        total : (bool), optional
            plot total concentrations. The default is True.
        ylabel_fontsize : (int), optional
            Size of the label for y-axis. The default is 15.
        **kwargs : (dict)
            Kwargs to pass to df.plot().

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        if conc_range:
            if not ylim:
                ylim = conc_range
            series = self.limit_conc_range(conc_range,reset_df=False)[1]
            df = self.limit_conc_range(conc_range,reset_df=False)[0]
        else:
            series = self.series
            df = self.df
        if not total:
            ax = df.plot(x='name',y='conc',kind='bar',ylim=ylim,logy=True,grid=True,
                          xlabel='Name , Charge',legend=None,**kwargs)
        else:
            ax = series.plot(kind='bar',logy=True,ylim=ylim,ylabel='Concentrations(cm$^{-3}$)',grid=True,**kwargs)
        ax.set_ylabel('Concentrations(cm$^{-3}$)',fontdict={'fontsize':ylabel_fontsize})
        return ax


    def limit_conc_range(self,conc_range=(1e10,1e40),reset_df=False):
        """
        Limit df to a concentration range

        Parameters
        ----------
        conc_range : (tuple), optional
            Range of concentrations to include in df. The default is (1e10,1e40).
        reset_df : (bool), optional
            Reset df attribute or return a new df. The default is False.

        Returns
        -------
        df, series : 
            DataFrame object, Series object.
        """
        if reset_df:
            self.df = self.df[self.df.conc.between(conc_range[0],conc_range[1])]
            self.series = self.series.between(conc_range[0],conc_range[1])
            return
        else:
            df = self.df.copy()
            series = self.series.copy()
            df = df[df.conc.between(conc_range[0],conc_range[1])]
            series = series[series.between(conc_range[0],conc_range[1])]
            return df , series
    
 
    
class ThermodynamicsPlotter:
    """
    Class to plot defect Thermodynamics data.
    """
    
    def __init__(self,xlim=(1e-20,1e08)):
        
        self.xlim = xlim
        

    def plot_pO2_vs_concentrations(self,thermodata,output='total',size=(12,8),fontsize=22,
                            xlim=None,ylim=None,show_unstable=True,colors=None,**kwargs):
        """
        Plot defect and carrier concentrations in a range of oxygen partial pressure.

        Parameters
        ----------
        thermodata: (ThermoData)
            ThermoData object that contains the thermodynamic data:
                partial_pressures : (list)
                    list with partial pressure values.
                defect_concentrations : (list or dict)
                    Defect concentrations in the same format as the output of DefectsAnalysis. 
                carrier_concentrations : (list)
                    List of tuples with intrinsic carriers concentrations (holes,electrons).
        output : (str), optional
            Type of output for defect concentrations:
                "all": The output is the concentration of every defect entry.
                "stable": The output is the concentration of the stable charge for every defect at each fermi level point.
                "total": The output is the sum of the concentration in every charge for each specie.
                The default is 'total'.
        size : (tuple), optional
            Size of the matplotlib figure. The default is (12,8).
        fontsize : (float), optional
            Size of font for matplotlib rcParams. The default is 22.
        xlim : (tuple), optional
            Range of x-axis. The default is (1e-20,1e08).
        ylim : (tuple), optional
            Range of y-axis. The default is None.
        show_unstable : (bool), optional
            Show regions where the system is unstable (at least one formation energy is negative).
        colors : (list)
            List of colors to use for plotting with matplotlib. If None the defaults are used.
        kwargs : 
            Kwargs to pass to DefectConcentrations.filter_concentrations(**kwargs).
            If provided, only the filtered concentrations will be plotted. If output
            is set to "total", only the filtered concentrations will be used to 
            compute the total concentration.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        td = thermodata
        p,dc,cc = td.partial_pressures,td.defect_concentrations,td.carrier_concentrations
        xlabel = 'Oxygen partial pressure (atm)'
        plt = self.plot_x_vs_concentrations(x=p,xlabel=xlabel,defect_concentrations=dc,
                                                 carrier_concentrations=cc,output=output,
                                                 size=size,fontsize=fontsize,xlim=xlim,ylim=ylim,
                                                 show_unstable=show_unstable,colors=colors,**kwargs)
        return plt
    

    def plot_pO2_vs_conductivity(self,partial_pressures,conductivities,new_figure=True,
                                 label=None,size=(12,8),xlim=None,ylim=None):
        """
        Plot conductivity as a function of the oxygen partial pressure.

        Parameters
        ----------
        partial_pressures : (list)
            list with partial pressure values.
        conductivities : (dict or list)
            If is a dict multiples lines will be plotted, with labels as keys and conductivity list
            as values. If is a list only one line is plotted with label taken from the "label" argument.
        new_figure : (bool), optional
            Initialize a new matplotlib figure. The default is True.
        label : (str), optional
            Label for the data. The default is None.
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
        xlabel = 'Oxygen partial pressure (atm)'
        plt = self.plot_x_vs_conductivity(x=partial_pressures,xlabel=xlabel,
                                          conductivities=conductivities,
                                          new_figure=new_figure,label=label,size=size,
                                          xlim=xlim,ylim=ylim)
        return plt    

    
    
    def plot_pO2_vs_fermi_level(self,partial_pressures,fermi_levels,band_gap,new_figure=True,
                         label=None,size=(12,8),xlim=None,ylim=None,colors=None):
        """
        Plot Fermi level as a function of the oxygen partial pressure.

        Parameters
        ----------
        partial_pressures : (list)
            list with partial pressure values.
        fermi_levels : (dict or list)
            If is a dict multiples lines will be plotted, with labels as keys and fermi level list
            as values. If is a list only one line is plotted with label taken from the "label" argument.
        band_gap : (float)
            Band gap of the bulk material.
        new_figure : (bool), optional
            Initialize a new matplotlib figure. The default is True.
        label : (str), optional
            Label for the data. The default is None.
        size : (tuple), optional
            Size of the matplotlib figure. The default is (12,8).
        xlim : (tuple), optional
            Range of x-axis. The default is (1e-20,1e08).
        ylim : (tuple), optional
            Range of y-axis. The default is None.
        colors : (list), optional
            List with colors for Fermi level data.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        xlabel = 'Oxygen partial pressure (atm)'
        plt = self.plot_x_vs_fermi_level(x=partial_pressures,xlabel=xlabel,
                                          fermi_levels=fermi_levels,band_gap=band_gap,
                                          new_figure=new_figure,label=label,size=size,
                                          xlim=xlim,ylim=ylim,colors=colors)
        return plt


    def plot_variable_species_vs_concentrations(self,thermodata,output='total',size=(12,8),fontsize=22,
                            xlim=None,ylim=None,show_unstable=True,colors=None,**kwargs):
        """
        Plot defect and carrier concentrations in a range of oxygen partial pressure.

        Parameters
        ----------
        thermodata: (ThermoData)
            ThermoData object that contains the thermodynamic data:
                variable_defect_specie : (str)
                    Name of variable defect species.
                variable_concentrations : (list)
                    List of concentrations of variable species. 
                defect_concentrations : (list or dict)
                    Defect concentrations in the same format as the output of DefectsAnalysis. 
                carrier_concentrations : (list)
                    List of tuples with intrinsic carriers concentrations (holes,electrons).
        output : (str), optional
            Type of output for defect concentrations:
                "all": The output is the concentration of every defect entry.
                "stable": The output is the concentration of the stable charge for every defect at each fermi level point.
                "total": The output is the sum of the concentration in every charge for each specie.
                The default is 'total'.
        size : (tuple), optional
            Size of the matplotlib figure. The default is (12,8).
        fontsize : (float), optional
            Size of font for matplotlib rcParams. The default is 22.
        xlim : (tuple), optional
            Range of x-axis. The default is (1e-20,1e08).
        ylim : (tuple), optional
            Range of y-axis. The default is None.
        show_unstable : (bool), optional
            Show regions where the system is unstable (at least one formation energy is negative).
        colors : (list)
            List of colors to use for plotting with matplotlib. If None the defaults are used.
        kwargs : 
            Kwargs to pass to DefectConcentrations.filter_concentrations(**kwargs).
            If provided, only the filtered concentrations will be plotted. If output
            is set to "total", only the filtered concentrations will be used to 
            compute the total concentration.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        td = thermodata
        c,dc,cc = td.variable_concentrations,td.defect_concentrations,td.carrier_concentrations
        dname = self._get_variable_defect_specie_label(td.variable_defect_specie)
        xlabel = '[%s] (cm$^{-3})$' %dname
        
        plt = self.plot_x_vs_concentrations(x=c,xlabel=xlabel,defect_concentrations=dc,
                                                 carrier_concentrations=cc,output=output,
                                                 size=size,fontsize=fontsize,xlim=xlim,ylim=ylim,
                                                 show_unstable=show_unstable,colors=colors,**kwargs)
        return plt


    def plot_variable_species_vs_conductivity(self,xlabel,variable_concentrations,conductivities,
                                              new_figure=True,label=None,size=(12,8),xlim=None,ylim=None):
        """
        Plot conductivity as a function of the oxygen partial pressure.

        Parameters
        ----------
        xlabel : (str)
            Label for concentration axis (cm^-3 is added).
        variable_concentrations : (list)
            List of concentrations of variable species. 
        conductivities : (dict or list)
            If is a dict multiples lines will be plotted, with labels as keys and conductivity list
            as values. If is a list only one line is plotted with label taken from the "label" argument.
        new_figure : (bool), optional
            Initialize a new matplotlib figure. The default is True.
        label : (str), optional
            Label for the data. The default is None.
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
        xlabel = '%s (cm$^{-3})$' %xlabel
        plt = self.plot_x_vs_conductivity(x=variable_concentrations,xlabel=xlabel,
                                          conductivities=conductivities,
                                          new_figure=new_figure,label=label,size=size,
                                          xlim=xlim,ylim=ylim)
        return plt
        

    def plot_variable_species_vs_fermi_level(self,xlabel,variable_concentrations,fermi_levels,band_gap,
                         new_figure=True,label=None,size=(12,8),xlim=None,ylim=None,colors=None):
        """
        Plot Fermi level as a function of the oxygen partial pressure.

        Parameters
        ----------
        xlabel : (str)
            Label for concentration axis (cm^-3 is added).
        variable_concentrations : (list)
            List of concentrations of variable species. 
        fermi_levels : (dict or list)
            If is a dict multiples lines will be plotted, with labels as keys and fermi level list
            as values. If is a list only one line is plotted with label taken from the "label" argument.
        band_gap : (float)
            Band gap of the bulk material.
        new_figure : (bool), optional
            Initialize a new matplotlib figure. The default is True.
        label : (str), optional
            Label for the data. The default is None.
        size : (tuple), optional
            Size of the matplotlib figure. The default is (12,8).
        xlim : (tuple), optional
            Range of x-axis. The default is (1e-20,1e08).
        ylim : (tuple), optional
            Range of y-axis. The default is None.
        colors : (list), optional
            List with colors for Fermi level data.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        xlabel = '%s (cm$^{-3})$' %xlabel
        plt = self.plot_x_vs_fermi_level(x=variable_concentrations,xlabel=xlabel,
                                          fermi_levels=fermi_levels,band_gap=band_gap,
                                          new_figure=new_figure,label=label,size=size,
                                          xlim=xlim,ylim=ylim,colors=colors)
        return plt



    def plot_x_vs_concentrations(self,x,xlabel,defect_concentrations,carrier_concentrations,
                                     output='total',size=(12,8),fontsize=22,
                                     xlim=None,ylim=None,show_unstable=True,colors=None,**kwargs):
        """
        Plot defect concentrations as a function of generic data on the x-axis.        
        
        Parameters
        ----------
        x : (list)
            List with data on x-axis.
        xlabel : (str)
            Label for x-axis.
        defect_concentrations : (list)
            List of DefectConcentrations objects.
        carrier_concentrations : (list)
            List of tuples with carrier concentrations (holes,electrons).
         output : (str), optional
             Type of output for defect concentrations:
                 "all": The output is the concentration of every defect entry.
                 "stable": The output is the concentration of the stable charge for every defect at each fermi level point.
                 "total": The output is the sum of the concentration in every charge for each specie.
                 The default is 'total'.
         size : (tuple), optional
             Size of the matplotlib figure. The default is (12,8).
         fontsize : (float), optional
             Size of font for matplotlib rcParams. The default is 22.
         xlim : (tuple), optional
             Range of x-axis. The default is (1e-20,1e08).
         ylim : (tuple), optional
             Range of y-axis. The default is None.
         show_unstable : (bool), optional
             Show regions where the system is unstable (at least one formation energy is negative).
         colors : (list)
             List of colors to use for plotting with matplotlib. If None the defaults are used.
         kwargs : 
             Kwargs to pass to DefectConcentrations.filter_concentrations(**kwargs).
             If provided, only the filtered concentrations will be plotted. If output
             is set to "total", only the filtered concentrations will be used to 
             compute the total concentration.

         Returns
         -------
         plt : 
             Matplotlib object.
         """        
        matplotlib.rcParams.update({'font.size': 22})
        if output == 'all' or output == 'stable':
            plt = self._plot_conc(x=x,defect_concentrations=defect_concentrations,
                             carrier_concentrations=carrier_concentrations,
                             output=output,size=size,colors=colors,**kwargs)
        elif output == 'total':
            plt = self._plot_conc_total(x=x,defect_concentrations=defect_concentrations,
                                   carrier_concentrations=carrier_concentrations,
                                   size=size,colors=colors,**kwargs)
        else:
            raise ValueError('The options for plot output are "all", "stable" or "total".')
            
        plt.xscale('log')
        plt.yscale('log')
        if show_unstable:
            stable = self._get_unstable_bool(defect_concentrations)
            ax = plt.gca()
            plt.fill_between(x, 0, 1, where=stable, alpha=0.3, transform=ax.get_xaxis_transform(), color='red')
        xlim = xlim if xlim else self.xlim
        plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)

        plt.xlabel(xlabel)
        plt.ylabel('Concentrations (cm$^{-3})$')
        plt.legend()
        plt.grid()
        
        return plt   

        

    def plot_x_vs_conductivity(self,x,xlabel,conductivities,new_figure=True,
                               label=None,size=(12,8),xlim=None,ylim=None):
        """
        Plot conductivity as a function of the oxygen partial pressure.

        Parameters
        ----------
        x : (list)
            List with data on x-axis.
        xlabel : (str)
            Label for x-axis.
        conductivities : (dict or list)
            If is a dict multiples lines will be plotted, with labels as keys and conductivity list
            as values. If is a list only one line is plotted with label taken from the "label" argument.
        new_figure : (bool), optional
            Initialize a new matplotlib figure. The default is True.
        label : (str), optional
            Label for the data. The default is None.
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
        matplotlib.rcParams.update({'font.size': 22})
        if new_figure:
            plt.figure(figsize=(size))
        if isinstance(conductivities,dict):
            for name,sigma in conductivities.items():
                plt.plot(x,sigma,linewidth=4,marker='s',label=name)
        else:
            sigma = conductivities
            plt.plot(x,sigma,linewidth=4,marker='s',label=label)
        plt.xscale('log')
        plt.yscale('log')
        xlim = xlim if xlim else self.xlim
        plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        plt.xlabel(xlabel)
        plt.ylabel('Conductivity (S/m)')
        plt.legend()
        if new_figure:
            plt.grid()
        
        return plt


    def plot_x_vs_fermi_level(self,x,xlabel,fermi_levels,band_gap,new_figure=True,
                         label=None,size=(12,8),xlim=None,ylim=None,colors=None):
        """
        Parameters
        ----------
        x : (list)
            List with data on x-axis.
        xlabel : (str)
            Label for x-axis.
        fermi_levels : (dict or list)
            If is a dict multiples lines will be plotted, with labels as keys and fermi level list
            as values. If is a list only one line is plotted with label taken from the "label" argument.
        band_gap : (float)
            Band gap of the bulk material.
        new_figure : (bool), optional
            Initialize a new matplotlib figure. The default is True.
        label : (str), optional
            Label for the data. The default is None.
        size : (tuple), optional
            Size of the matplotlib figure. The default is (12,8).
        xlim : (tuple), optional
            Range of x-axis. The default is (1e-20,1e08).
        ylim : (tuple), optional
            Range of y-axis. The default is None.
        colors : (list), optional
            List with colors for Fermi level data.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        ylim = ylim if ylim else (-0.5, band_gap+0.5)
        matplotlib.rcParams.update({'font.size': 22})
        if new_figure:
            plt.figure(figsize=(size))
        if isinstance(fermi_levels,dict):
            for name,mue in fermi_levels.items():
                clr = colors[list(fermi_levels.keys()).index(name)] if colors else None
                plt.plot(x,mue,linewidth=4,label=name,color=clr)
        else:
            mue = fermi_levels
            clr = colors[0] if colors else None
            plt.plot(x,mue,linewidth=4,label=label,color=clr)
        plt.xscale('log')
        xlim = xlim if xlim else self.xlim
        plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        plt.hlines(0, xlim[0], xlim[1],color='k')
        plt.text(xlim[1]+(xlim[1]-xlim[0]),-0.05,'VB')
        plt.text(xlim[1]+(xlim[1]-xlim[0]),band_gap-0.05,'CB')
        plt.hlines(band_gap, xlim[0], xlim[1],color='k')
        plt.axhspan(ylim[0], 0, alpha=0.2,color='k')
        plt.axhspan(band_gap, ylim[1], alpha=0.2,color='k')

        plt.xlabel(xlabel)
        plt.ylabel('Electron chemical potential (eV)')
        plt.legend()
        if new_figure:
            plt.grid()
        
        return plt
    
    
    def _get_variable_defect_specie_label(self,variable_defect_specie):
        try:
            dname = get_defect_name_from_string(variable_defect_specie)
            return dname.symbol
        except:
            return variable_defect_specie

    
    
    def _plot_conc(self,x,defect_concentrations,carrier_concentrations,output,size,colors,**kwargs):
        
        plt.figure(figsize=size)
        dc = defect_concentrations if output != 'stable' else [c.stable for c in defect_concentrations] 
        if kwargs:
            dc = [c.filter_concentrations(**kwargs) for c in dc] #filter concentrations based on kwargs
        h = [cr[0] for cr in carrier_concentrations] 
        n = [cr[1] for cr in carrier_concentrations]
        previous_charge = None
        for i in range(0,len(dc[0])):
                conc = [c[i].conc for c in dc]
                charges = [c[i].charge for c in dc]
                label_txt = dc[0][i].name.symbol
                if output == 'all':
                    label_txt = format_legend_with_charge_number(label_txt,dc[0][i].charge)
                elif output == 'stable':
                    for q in charges:
                        if q != previous_charge:
                            previous_charge = q
                            label_charge = '+' + str(int(q)) if q > 0 else str(int(q))
                            index = charges.index(q)
                            plt.text(x[index],conc[index],label_charge,clip_on=True)
                color = colors[i] if colors else None
                plt.plot(x,conc,label=label_txt,linewidth=4,color=color)
        plt.plot(x,h,label='$n_{h}$',linestyle='--',color='r',linewidth=4)
        plt.plot(x,n,label='$n_{e}$',linestyle='--',color='b',linewidth=4)
        return plt
    
    
    def _plot_conc_total(self,x,defect_concentrations,carrier_concentrations,size,colors,**kwargs):
        
        dc = defect_concentrations
        if kwargs:
            dc = [c.filter_concentrations(**kwargs) for c in dc]
        plt.figure(figsize=size)
        h = [cr[0] for cr in carrier_concentrations] 
        n = [cr[1] for cr in carrier_concentrations]
        for name in dc[0].names:
            conc = [c.total[name] for c in dc]
            label_txt = name.symbol
            color = colors[dc[0].names.index(name)] if colors else None
            plt.plot(x,conc,label=label_txt,linewidth=4,color=color)
        plt.plot(x,h,label='$n_{h}$',linestyle='--',color='r',linewidth=4)
        plt.plot(x,n,label='$n_{e}$',linestyle='--',color='b',linewidth=4)
        return plt
    
    
    def _get_unstable_bool(self,defect_concentrations):
        slist = []
        for dc in defect_concentrations:
            unstable = False
            for d in dc:
                if 'stable' in vars(d).keys() and d.stable == False:
                    unstable = True
            slist.append(unstable)
        
        return slist
        
    
    
    
