#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 15:24:29 2021

@author: villa
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from .defects import Defect, format_legend_with_charge_number, get_defect_from_string



def plot_formation_energies(entries,
                            chemical_potentials,
                            vbm,
                            band_gap,
                            temperature=0,
                            xlim=None,
                            ylim=None,
                            title=None,
                            fermi_level=None,
                            grid=True,
                            figsize=(6,6),
                            fontsize=12,
                            show_legend=True,
                            format_legend=True,
                            get_subplot=False,
                            subplot_settings=None,
                            **eform_kwargs):
    """
    Produce defect Formation energy vs Fermi energy plot.

    Parameters
    ----------
    entries : (list) 
        List of entries to calculate.
    chemical_potentials : (dict)
        Dictionary with chemical potentials of the elements {'element':chempot}.
    vbm : (float)
        Valence band maximum in eV.
    band_gap : (float)
        Band gap of bulk material in eV.
    temperature : (float)
        Temperature in K. If no custom formation energy is provided, this arg has no effect. 
    xlim : (tuple)
        Tuple (min,max) giving the range of the x (fermi energy) axis.
    ylim : (tuple)
        Tuple (min,max) giving the range for the formation energy axis.
    title : (str)
        Title of the figure.
    fermi_level : (float)
        Plot Fermi energy position with a vertical line.
    grid : (bool)
        Show grid.
    figsize : (float or tuple)
        Figure size.
    fontsize : (float)
        Font size.
    show_legend  : (Bool)
        Show legend.
    format_legend : (bool)
        Get latex-like legend based on the name of defect entries.
    get_subplot : (bool)
        Get subplot.
    subplot_settings:
        List with integers for subplot setting on matplotlib (plt.subplot(nrows,ncolumns,index)). 
    eform_kwargs : (dict)
        Kwargs to pass to `entry.formation_energy`.

    Returns
    -------
        matplotlib object
    """
    from .analysis import DefectsAnalysis
    import matplotlib.pyplot as plt
    
    matplotlib.rcParams.update({'font.size': fontsize}) 
    formation_energies = DefectsAnalysis(entries,vbm,band_gap).formation_energies(chemical_potentials=chemical_potentials,
                                                                                    fermi_level=0,
                                                                                    temperature=temperature,
                                                                                    entries=entries,
                                                                                    **eform_kwargs)
    if xlim == None:
        xlim = (-0.5,band_gap+0.5)        
    npoints = 200
    step = abs(xlim[1]+0.1-xlim[0])/npoints
    x = np.arange(xlim[0],xlim[1]+0.1,step)
    
    if get_subplot:
        if subplot_settings[2] == 1:
            plt.figure(figsize=figsize)               
        plt.subplot(subplot_settings[0],subplot_settings[1],subplot_settings[2])
    else:
        plt.figure(figsize=figsize)
        
    for name in formation_energies:
        energy = np.zeros(len(x))
        emin = np.zeros(len(x))
        x_star = []
        y_star = []
        q_previous = None
        for i in range(0,npoints):
            emin[i] = 1e40
            for d in formation_energies[name]:
                q = d[0]
                e0 = d[1]
                energy[i] = e0 + q*x[i]
                # finding most stable charge state
                if energy[i] < emin[i]:
                    emin[i] = energy[i]
                    q_stable = q
            # getting data to plot transition levels        
            if q_stable != q_previous:
                if q_previous != None:
                    x_star.append(x[i])
                    y_star.append(emin[i])
                q_previous = q_stable       

        if format_legend:
            label_txt = get_defect_from_string(name).symbol
        else:
            label_txt = name            

        plt.plot(x,emin,label=label_txt,linewidth=3)
        plt.scatter(x_star,y_star,s=120,marker='*')
                    
    plt.axvline(x=0.0, linestyle='-', color='k', linewidth=2)  # black dashed lines for gap edges
    plt.axvline(x=band_gap, linestyle='-', color='k',
                linewidth=2)        
    if fermi_level:
        plt.axvline(x=fermi_level, linestyle='dashed', color='k', linewidth=1.5, label='$\mu _{e}$')                
    # shaded areas
    plt.axvspan(xlim[0], 0, facecolor='k', alpha=0.2)
    plt.axvspan(band_gap, xlim[1]+0.1, facecolor='k', alpha=0.2)
    plt.hlines(0,xlim[0],xlim[1]+0.1,colors='k',linestyles='dashed',alpha=0.5)
    plt.xlim(xlim)
    if ylim: 
        plt.ylim(ylim) 
    plt.xlabel('Fermi level (eV)')
    plt.ylabel('Formation energy (eV)')
    if title:
        plt.title(title)
    if show_legend:    
        plt.legend()
    if grid:
        plt.grid()

    return plt


def plot_binding_energies(entries,
                          vbm,
                          band_gap,
                          temperature=0,
                          names=None,
                          xlim=None,
                          ylim=None,
                          figsize=(6,6),
                          fontsize=18,
                          format_legend=True,
                          **eform_kwargs):
    """
    Plot binding energies for complex of defects as a function of the fermi level

    Parameters
    ----------
    entries : (list) 
        List of entries to calculate.
    vbm : (float)
        Valence band maximum in eV.
    band_gap : (float)
        Band gap of bulk material in eV.
    temperature : (float)
        Temperature in K. If no custom formation energy is provided, this arg has no effect.
    names : (list)
        List of strings with names of DefectEntry. If None, all defect complexes are plotted.
    xlim : (tuple)
        Tuple (min,max) giving the range of the x (fermi energy) axis
    ylim : (tuple)
        Tuple (min,max) giving the range for the formation energy axis
    figsize : (tuple)
        Figure size.
    fontsize : (float)
        Font size.
    format_legend : (bool)
        Bool for getting latex-like legend based on the name of defect entries.
    eform_kwargs : (dict)
        Kwargs to pass to `entry.formation_energy`.

    Returns
    -------
        matplotlib object
    """         
    from .analysis import DefectsAnalysis
    
    da = DefectsAnalysis(entries, vbm, band_gap)
    plt.figure(figsize=figsize)
    matplotlib.rcParams.update({'font.size': fontsize}) 
    if xlim==None:
        xlim = (-0.5,da.band_gap+0.5)
    # building array for x values (fermi level)    
    ef = np.arange(xlim[0],xlim[1]+0.1,(xlim[1]-xlim[0])/200)        
    binding_energy = np.zeros(len(ef))        
    if not names:
        names = []
        for e in da.entries:
            if e.defect_type == 'DefectComplex':
                if e.name not in names:
                    names.append(e.name)   
    if not names:
        raise ValueError('No DefectComplex entries found')  

    # getting binding energy at different fermi levels for every name in list
    for name in names:
        label = da.select_entries(names=[name])[0].symbol if format_legend else name
        for i in range(0,len(ef)):
            binding_energy[i] = da.binding_energy(
                                        name=name,
                                        fermi_level=ef[i],
                                        temperature=temperature,
                                        **eform_kwargs)
                        
        plt.plot(ef,binding_energy, linewidth=2.5*(figsize[1]/figsize[0]),label=label)
        
    plt.axvline(x=0.0, linestyle='-', color='k', linewidth=2)  # black lines for gap edges
    plt.axvline(x=da.band_gap, linestyle='-', color='k',
                linewidth=2)        
    # shaded areas
    plt.axvspan(xlim[0], 0, facecolor='k', alpha=0.2)
    plt.axvspan(da.band_gap, xlim[1], facecolor='k', alpha=0.2)
    plt.hlines(0,xlim[0],xlim[1],colors='k',linestyles='dashed',alpha=0.5)
    plt.legend()
    plt.xlim(xlim)
    if ylim: 
        plt.ylim(ylim) 
    plt.xlabel('Fermi level (eV)')
    plt.ylabel('Binding energy (eV)')
    
    return plt


def plot_charge_transition_levels(entries,
                                  vbm,
                                  band_gap,
                                  temperature=0,
                                  ylim=None,
                                  figsize=(8,8),
                                  fontsize=16,
                                  fermi_level=None,
                                  format_legend=True,
                                  get_integers=True,
                                  **eform_kwargs):
    """
    Plotter for the charge transition levels.

    Parameters
    ----------
    entries : (list) 
        List of entries to calculate.
    vbm : (float)
        Valence band maximum in eV.
    band_gap : (float)
        Band gap of bulk material in eV.
    temperature : (float)
        Temperature in K. If no custom formation energy is provided, this arg has no effect.
    ylim : (tuple)
        y-axis limits.
    figsize : (tuple)
        Figure size.
    fontsize : (float)
        Font size.
    fermi_level : (float)
        Plot Fermi energy position.
    format_legend : (bool)
        Bool for getting latex-like legend based on the name of defect entries.
    get_integers : (bool)
        Get charge transition levels as integers.
    eform_kwargs : (dict)
        Kwargs to pass to `entry.formation_energy`.

    Returns
    -------
        matplotlib object    
    """        
    from .analysis import DefectsAnalysis
    
    da = DefectsAnalysis(entries, vbm, band_gap)
    plt.figure(figsize=figsize)         
    if ylim == None:
        ylim = (-0.5,da.band_gap +0.5)        
    charge_transition_levels = da.charge_transition_levels(
                                                temperature=temperature,
                                                entries=entries,
                                                get_integers=get_integers,
                                                **eform_kwargs)
    
    number_defects = len(charge_transition_levels)   
    x_max = 10
    interval = x_max/(number_defects + 1)
    x = np.arange(0,x_max,interval)
    # position of x labels
    x_ticks_positions = []
    for i in range(0,len(x)-1):
        x_ticks_positions.append((x[i+1]-x[i])/2 + x[i])            
    x_ticks_labels = []
    for name in charge_transition_levels:
        x_ticks_labels.append(name)        
    # draw vertical lines to separte defect types
    for i in x:
        plt.axvline(x=i, linestyle='-', color='k', linewidth=1.2, alpha=1, zorder=1)
    xlim = (x[0],x[-1])        
    #VBM and CBM shaded
    plt.axhspan(ylim[0], 0, facecolor='grey', alpha=0.9, zorder=2)
    plt.axhspan(da.band_gap,ylim[1], facecolor = 'grey', alpha=0.9, zorder=2)                
    # plot CTL
    for i in range(0,len(x_ticks_labels)):
        name = x_ticks_labels[i]
        for ctl in charge_transition_levels[name]:
            energy = ctl[2]
            plt.hlines(energy,x[i],x[i+1],colors='k',linewidth=2.25, zorder=3)
            charge1 = '+' + str(ctl[1]) if ctl[1] > 0 else str(ctl[1])
            charge2 = '+' + str(ctl[0]) if ctl[0] > 0 else str(ctl[0])
            label_charge = '(' + charge2 + '/' + charge1 + ')'
            font_space = abs(ylim[1]-ylim[0]) / 100
            if energy < ylim[1] and energy > ylim[0]:
                plt.text(x[i]+(interval/2)*2/number_defects ,energy+font_space,label_charge,fontsize=fontsize)        
    # format latex-like legend
    if format_legend:    
         for name in x_ticks_labels:            
            x_ticks_labels[x_ticks_labels.index(name)] = get_defect_from_string(name).symbol               
    if fermi_level:
        plt.axhline(y=fermi_level, linestyle='dashed', color='k', linewidth=1.5, label='$\mu _{e}$')   
    
    plt.text(x[-1]+interval/8,-0.3,'VB',fontsize=25*(fontsize/16))
    plt.text(x[-1]+interval/8,da.band_gap+0.2,'CB',fontsize=25*(fontsize/16))
    plt.xticks(ticks=x_ticks_positions,labels=x_ticks_labels,fontsize = (25-number_defects)*(fontsize/16))
    plt.tick_params(axis='x',length=0,width=0)
    plt.yticks(fontsize=16*(fontsize/16))
    plt.xlim(xlim)  
    plt.ylim(ylim)
    plt.ylabel('Energy(eV)',fontsize=20*(fontsize/16))  
    
    return plt    


def plot_pO2_vs_concentrations(
                            thermodata,
                            output='total',
                            figsize=(8,8),
                            fontsize=22,
                            xlim=(1e-20,1e10),
                            ylim=None,
                            show_unstable=True,
                            colors=None,
                            **kwargs):
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
    output : (str)
        Type of output for defect concentrations:
            "all": The output is the concentration of every defect entry.
            "stable": The output is the concentration of the stable charge for every defect at each fermi level point.
            "total": The output is the sum of the concentration in every charge for each specie.
            The default is 'total'.
    figsize : (tuple)
        Size of the matplotlib figure.
    fontsize : (float)
        Size of font for matplotlib rcParams.
    xlim : (tuple)
        Range of x-axis.
    ylim : (tuple)
        Range of y-axis. 
    show_unstable : (bool)
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
    plt = plot_x_vs_concentrations(
                                x=p,
                                xlabel=xlabel,
                                defect_concentrations=dc,
                                carrier_concentrations=cc,
                                output=output,
                                figsize=figsize,
                                fontsize=fontsize,
                                xlim=xlim,
                                ylim=ylim,
                                show_unstable=show_unstable,
                                colors=colors,
                                **kwargs)
    return plt


def plot_pO2_vs_conductivity(
                            partial_pressures,
                            conductivities,
                            new_figure=True,
                            label=None,
                            figsize=(8,8),
                            xlim=(1e-20,1e10),
                            ylim=None):
    """
    Plot conductivity as a function of the oxygen partial pressure.

    Parameters
    ----------
    partial_pressures : (list)
        list with partial pressure values.
    conductivities : (dict or list)
        If is a dict multiples lines will be plotted, with labels as keys and conductivity list
        as values. If is a list only one line is plotted with label taken from the "label" argument.
    new_figure : (bool)
        Initialize a new matplotlib figure.
    label : (str)
        Label for the data.
    figsize : (tuple)
        Size of the matplotlib figure.
    xlim : (tuple)
        Range of x-axis.
    ylim : (tuple)
        Range of y-axis.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    xlabel = 'Oxygen partial pressure (atm)'
    plt = plot_x_vs_conductivity(
                                x=partial_pressures,
                                xlabel=xlabel,
                                conductivities=conductivities,
                                new_figure=new_figure,
                                label=label,
                                figsize=figsize,
                                xlim=xlim,
                                ylim=ylim)
    return plt    



def plot_pO2_vs_fermi_level(
                            partial_pressures,
                            fermi_levels,
                            band_gap,
                            new_figure=True,
                            label=None,
                            figsize=(8,8),
                            xlim=(1e-20,1e10),
                            ylim=None,
                            colors=None):
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
    new_figure : (bool)
        Initialize a new matplotlib figure.
    label : (str)
        Label for the data.
    figsize : (tuple)
        Size of the matplotlib figure.
    xlim : (tuple)
        Range of x-axis.
    ylim : (tuple)
        Range of y-axis.
    colors : (list)
        List with colors for Fermi level data.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    xlabel = 'Oxygen partial pressure (atm)'
    plt = plot_x_vs_fermi_level(
                                x=partial_pressures,
                                xlabel=xlabel,
                                fermi_levels=fermi_levels,
                                band_gap=band_gap,
                                new_figure=new_figure,
                                label=label,
                                figsize=figsize,
                                xlim=xlim,
                                ylim=ylim,
                                colors=colors)
    return plt


def plot_variable_species_vs_concentrations(
                                        thermodata,
                                        output='total',
                                        figsize=(8,8),
                                        fontsize=22,
                                        xlim=(1e-20,1e10),
                                        ylim=None,
                                        show_unstable=True,
                                        colors=None,
                                        **kwargs):
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
    output : (str)
        Type of output for defect concentrations:
            "all": The output is the concentration of every defect entry.
            "stable": The output is the concentration of the stable charge for every defect at each fermi level point.
            "total": The output is the sum of the concentration in every charge for each specie.
            The default is 'total'.
    figsize : (tuple)
        Size of the matplotlib figure.
    fontsize : (float)
        Size of font for matplotlib rcParams.
    xlim : (tuple)
        Range of x-axis.
    ylim : (tuple)
        Range of y-axis.
    show_unstable : (bool)
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
    dname = _get_variable_defect_specie_label(td.variable_defect_specie)
    xlabel = '[%s] (cm$^{-3})$' %dname
    
    plt = plot_x_vs_concentrations(
                                x=c,
                                xlabel=xlabel,
                                defect_concentrations=dc,
                                carrier_concentrations=cc,
                                output=output,
                                figsize=figsize,
                                fontsize=fontsize,
                                xlim=xlim,
                                ylim=ylim,
                                show_unstable=show_unstable,
                                colors=colors,
                                **kwargs)
    return plt


def plot_variable_species_vs_conductivity(
                                        xlabel,
                                        variable_concentrations,
                                        conductivities,
                                        new_figure=True,
                                        label=None,
                                        figsize=(8,8),
                                        xlim=(1e-20,1e10),
                                        ylim=None):
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
    new_figure : (bool)
        Initialize a new matplotlib figure.
    label : (str)
        Label for the data.
    figsize : (tuple)
        Size of the matplotlib figure.
    xlim : (tuple)
        Range of x-axis.
    ylim : (tuple)
        Range of y-axis.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    xlabel = '%s (cm$^{-3})$' %xlabel
    plt = plot_x_vs_conductivity(
                                x=variable_concentrations,
                                xlabel=xlabel,
                                conductivities=conductivities,
                                new_figure=new_figure,
                                label=label,
                                figsize=figsize,
                                xlim=xlim,
                                ylim=ylim)
    return plt
    

def plot_variable_species_vs_fermi_level(
                                        xlabel,
                                        variable_concentrations,
                                        fermi_levels,
                                        band_gap,
                                        new_figure=True,
                                        label=None,
                                        figsize=(8,8),
                                        xlim=(1e-20,1e10),
                                        ylim=None,
                                        colors=None):
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
    new_figure : (bool)
        Initialize a new matplotlib figure.
    label : (str)
        Label for the data.
    figsize : (tuple)
        Size of the matplotlib figure.
    xlim : (tuple)
        Range of x-axis.
    ylim : (tuple)
        Range of y-axis.
    colors : (list)
        List with colors for Fermi level data.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    xlabel = '%s (cm$^{-3})$' %xlabel
    plt = plot_x_vs_fermi_level(
                                x=variable_concentrations,
                                xlabel=xlabel,
                                fermi_levels=fermi_levels,
                                band_gap=band_gap,
                                new_figure=new_figure,
                                label=label,
                                figsize=figsize,
                                xlim=xlim,
                                ylim=ylim,
                                colors=colors)
    return plt



def plot_x_vs_concentrations(
                            x,
                            xlabel,
                            defect_concentrations,
                            carrier_concentrations,
                            output='total',
                            figsize=(8,8),
                            fontsize=22,
                            xlim=(1e-20,1e10),
                            ylim=None,
                            show_unstable=True,
                            colors=None,
                            **kwargs):
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
        output : (str)
            Type of output for defect concentrations:
                "all": The output is the concentration of every defect entry.
                "stable": The output is the concentration of the stable charge for every defect at each fermi level point.
                "total": The output is the sum of the concentration in every charge for each specie.
                The default is 'total'.
        figsize : (tuple)
            Size of the matplotlib figure.
        fontsize : (float)
            Size of font for matplotlib rcParams. 
        xlim : (tuple)
            Range of x-axis. 
        ylim : (tuple)
            Range of y-axis.
        show_unstable : (bool)
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
        plt = _plot_conc(
                        x=x,
                        defect_concentrations=defect_concentrations,
                        carrier_concentrations=carrier_concentrations,
                        output=output,figsize=figsize,colors=colors,**kwargs)
        
    elif output == 'total':
        plt = _plot_conc_total(
                            x=x,
                            defect_concentrations=defect_concentrations,
                            carrier_concentrations=carrier_concentrations,
                            figsize=figsize,colors=colors,**kwargs)
    else:
        raise ValueError('The options for plot output are "all", "stable" or "total".')
        
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)

    plt.xlabel(xlabel)
    plt.ylabel('Concentrations (cm$^{-3})$')
    plt.legend()
    plt.grid()
    
    return plt   

    

def plot_x_vs_conductivity(
                        x,
                        xlabel,
                        conductivities,
                        new_figure=True,
                        label=None,
                        figsize=(8,8),
                        xlim=(1e-20,1e10),
                        ylim=None):
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
    new_figure : (bool)
        Initialize a new matplotlib figure. 
    label : (str)
        Label for the data.
    figsize : (tuple)
        Size of the matplotlib figure.
    xlim : (tuple)
        Range of x-axis.
    ylim : (tuple)
        Range of y-axis.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    matplotlib.rcParams.update({'font.size': 22})
    if new_figure:
        plt.figure(figsize=figsize)
    if isinstance(conductivities,dict):
        for name,sigma in conductivities.items():
            plt.plot(x,sigma,linewidth=4,marker='s',label=name)
    else:
        sigma = conductivities
        plt.plot(x,sigma,linewidth=4,marker='s',label=label)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.xlabel(xlabel)
    plt.ylabel('Conductivity (S/m)')
    plt.legend()
    if new_figure:
        plt.grid()
    
    return plt


def plot_x_vs_fermi_level(
                        x,
                        xlabel,
                        fermi_levels,
                        band_gap,
                        new_figure=True,
                        label=None,
                        figsize=(8,8),
                        xlim=(1e-20,1e10),
                        ylim=None,
                        colors=None):
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
    new_figure : (bool)
        Initialize a new matplotlib figure.
    label : (str)
        Label for the data.
    figsize : (tuple)
        Size of the matplotlib figure.
    xlim : (tuple)
        Range of x-axis. 
    ylim : (tuple)
        Range of y-axis.
    colors : (list),
        List with colors for Fermi level data.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    ylim = ylim if ylim else (-0.5, band_gap+0.5)
    matplotlib.rcParams.update({'font.size': 20})
    if new_figure:
        plt.figure(figsize=figsize)
    if isinstance(fermi_levels,dict):
        for name,mue in fermi_levels.items():
            clr = colors[list(fermi_levels.keys()).index(name)] if colors else None
            plt.plot(x,mue,linewidth=4,label=name,color=clr)
    else:
        mue = fermi_levels
        clr = colors[0] if colors else None
        plt.plot(x,mue,linewidth=4,label=label,color=clr)
    plt.xscale('log')
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
    if type(fermi_levels) == dict:
        plt.legend()
    if new_figure:
        plt.grid()
    
    return plt


def _get_variable_defect_specie_label(variable_defect_specie):
    try:
        defect = get_defect_from_string(variable_defect_specie)
        return defect.symbol
    except:
        return variable_defect_specie



def _plot_conc(
            x,
            defect_concentrations,
            carrier_concentrations,
            output,
            figsize,
            colors,
            **kwargs):
    
    plt.figure(figsize=figsize)
    dc = defect_concentrations if output != 'stable' else [c.stable for c in defect_concentrations] 
    if kwargs:
        dc = [c.filter_concentrations(**kwargs) for c in dc] #filter concentrations based on kwargs
    h = [cr[0] for cr in carrier_concentrations] 
    n = [cr[1] for cr in carrier_concentrations]
    previous_charge = None
    for i in range(0,len(dc[0])):
            conc = [c[i].conc for c in dc]
            charges = [c[i].charge for c in dc]
            label_txt = get_defect_from_string(dc[0][i].name).symbol
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


def _plot_conc_total(
                    x,
                    defect_concentrations,
                    carrier_concentrations,
                    figsize,
                    colors,
                    **kwargs):
    
    dc = defect_concentrations
    if kwargs:
        dc = [c.filter_concentrations(**kwargs) for c in dc]
    plt.figure(figsize=figsize)
    h = [cr[0] for cr in carrier_concentrations] 
    n = [cr[1] for cr in carrier_concentrations]
    for name in dc[0].names:
        conc = [c.total[name] for c in dc]
        label_txt = get_defect_from_string(name).symbol
        color = colors[dc[0].names.index(name)] if colors else None
        plt.plot(x,conc,label=label_txt,linewidth=4,color=color)
    plt.plot(x,h,label='$n_{h}$',linestyle='--',color='r',linewidth=4)
    plt.plot(x,n,label='$n_{e}$',linestyle='--',color='b',linewidth=4)
    return plt


def _get_unstable_bool(defect_concentrations):
    slist = []
    for dc in defect_concentrations:
        unstable = False
        for d in dc:
            if 'stable' in vars(d).keys() and d.stable == False:
                unstable = True
        slist.append(unstable)
    
    return slist
        
    

class DefectConcentrationsPlotter:

    def __init__(self,concentrations,format_names=True):
        """
        Class to visualize concentrations dictionaries with pd.DataFrame and 
        pd.Series (fortotal concentrations).

        Parameters
        ----------
        concentrations : (list or dict)
            Concentrations can be in different formats: list of dict for "all" and 
            "stable_charges" and dict for "total".
        format_names : (bool)
            Format names with latex symbols.
        """
        # to be fixed
        conc_dict = []
        for c in concentrations:
            d = {'charge':c.charge,'conc':c.conc}

            if format_names:
                name = get_defect_from_string(c.name).symbol
            else:
                name = c.name
            d['name'] = format_legend_with_charge_number(name,c.charge)
            conc_dict.append(d)
        
        conc_total_dict = {}
        for dn,conc in concentrations.total.items():
            if format_names:
                name = get_defect_from_string(dn).symbol
            else:
                name = dn
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
        conc_range : (tuple)
            Range of concentrations to include in df.
        ylim : (tuple)
            Limit of y-axis in plot. If None conc_range is used.
        total : (bool)
            plot total concentrations.
        ylabel_fontsize : (int)
            Size of the label for y-axis.
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
        conc_range : (tuple)
            Range of concentrations to include in df.
        reset_df : (bool)
            Reset df attribute or return a new df.

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
    
