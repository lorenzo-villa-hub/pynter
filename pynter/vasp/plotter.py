#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 11:19:30 2021

@author: villa
"""
from pymatgen.electronic_structure.plotter import BSPlotter,BSDOSPlotter,DosPlotter, BSPlotterProjected  #original pymatgen
#from pynter.vasp.pmg.pmg_electronic_structure_plotter import BSPlotter,BSDOSPlotter,DosPlotter, BSPlotterProjected #modified version 
from pymatgen.electronic_structure.bandstructure import BandStructure


def plot_bs(bs,**kwargs):
    """
    Get BS plot with pymatgen.

    Parameters
    ----------
    bs : 
        BandStructureSymmLine object, most likely generaten from Vasprun or BSVasprun.
    **kwargs : (dict)
        Arguments for the get_plot function in BSPlotter in pymatgen.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    plotter = BSPlotter(bs)
    plt = plotter.get_plot(**kwargs)
    
    return plt


def plot_bs_elt_projected(bs,legend=True,**kwargs):
    """
    Get element projected BS color plot with pymatgen. Note that vasprun.xml needs to
    be parsed with parse_projected_eigen set to True.

    Parameters
    ----------
    bs : 
        BandStructureSymmLine object, most likely generated from Vasprun or BSVasprun.
    legend: (bool or dict)
        Get colors legend. If is a dict it is passed as kwargs in plt.legend(**kwargs)
    **kwargs : (dict)
        Arguments for the get_elt_projected_plots_color function in BSPlotter in pymatgen.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    plotter = BSPlotterProjected(bs)
    plt = plotter.get_elt_projected_plots_color(**kwargs)
    
    import matplotlib.patches as mpatches
    colors = ['red','green','blue']
    patches = []
    if legend:
        if 'elt_ordered' not in kwargs.keys():
            elt_ordered = plotter._bs.structure.composition.elements
        for elt in elt_ordered:
            el = elt.symbol
            patches.append(mpatches.Patch(color=colors[elt_ordered.index(elt)],label=el))
    if type(legend) == dict:
        plt.legend(handles=patches,**legend)
    else:
        plt.legend(handles=patches)
    return plt


def plot_dos_bs(dos,bs,draw_fermi=True,**kwargs):
    """
    Plot DOS and BS with pymatgen's BSDOSPlotter class.

    Parameters
    ----------
    dos : 
        Dos object.
    bs : 
        BandStructureSymmLine object, most likely generated from Vasprun or BSVasprun..
    **kwargs : (dict)
        Arguments to pass to BSDOSPlotter class.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    if 'bs_projection' not in kwargs.keys(): # setting default projections to None
        kwargs['bs_projection'] = None
    if 'dos_projection' not in kwargs.keys():
        kwargs['dos_projection'] = None
    if draw_fermi:
        plt = BSDOSPlotter(**kwargs).get_plot(bs,dos) # ensure that also original pymatgen function works
    else:
        plt = BSDOSPlotter(**kwargs).get_plot(bs,dos,draw_fermi=draw_fermi)
        
        
    return plt


def plot_dos_dict(complete_dos,d,xlim=(-3,3),**kwargs):
    """
    Plot DOS dict with pymatgen. The dict is most likely generated with the methods of the
    pymatgen CompleteDos class.

    Parameters
    ----------
    complete_dos : 
        CompleteDos object.
    d : TYPE
        Dict of DOS.
    xlim : (tuple), optional
        Limits for x axis. The default is (-3,3).
    **kwargs : (dict)
        Arguments for DosPlotter in pymatgen.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    
    dos_plotter = DosPlotter(**kwargs)
    for k in d:
        dos_plotter.add_dos(k,d[k])
    eg = complete_dos.get_gap()
    plt = dos_plotter.get_plot(xlim=(xlim[0],eg+xlim[1]))
    
    return plt


def plot_element_dos(complete_dos,xlim=(-3,3),**kwargs):
    """
    Plot element resolved DOS with pymatgen

    Parameters
    ----------
    complete_dos : 
        CompleteDos object.
    xlim : (tuple), optional
        Limits for x axis. The default is (-3,3).
    **kwargs : (dict)
        Arguments for DosPlotter in pymatgen.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    d = complete_dos.get_element_dos()        
    plt = plot_dos_dict(complete_dos, d, xlim, **kwargs)
    
    return plt


def plot_spd_dos(complete_dos,xlim=(-3,3),**kwargs):
    """
    Plot spd resolved DOS with pymatgen

    Parameters
    ----------
    complete_dos : 
        CompleteDos object.
    xlim : (tuple), optional
        Limits for x axis. The default is (-3,3).
    **kwargs : (dict)
        Arguments for DosPlotter in pymatgen.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    d = complete_dos.get_spd_dos()        
    plt = plot_dos_dict(complete_dos, d, xlim, **kwargs)
    
    return plt


def plot_dos_bs_custom(bs,dos,ylim=(-4,7)):
    
    from pynter.vasp.pmg_electronic_structure_plotter import BSPlotter, DosPlotter
    # PLOT BS
    plt = BSPlotter(bs).get_plot(ylim=ylim,get_subplot=True)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
        
    # PLOT DOS
    partial_dos = dos.get_spd_dos()
    dos_plotter = DosPlotter()
    dos_plotter.add_dos('Total_dos',dos)
    for orbital in partial_dos:
        dos_plotter.add_dos(orbital,partial_dos[orbital])
    plt = dos_plotter.get_plot(xlim=ylim,get_subplot=True)    

    # modify pymatgen output to flip graph
    ax = plt.gca()
    ymax = 0
    fermi_lines = []
    for i in range(0,len(ax.lines)):
        x = ax.lines[i].get_ydata()
        y = ax.lines[i].get_xdata()
        for j in range(0,len(x)):
            if x[j] > ylim[0] and x[j] < ylim[1]:
                if x[j] > ymax:
                    ymax = x[j]
        ax.lines[i].set_xdata(x)
        ax.lines[i].set_ydata(y)
        if ax.lines[i].get_linestyle() == '--':
            fermi_lines.append(ax.lines[i])
            
    ax.lines = [l for l in ax.lines if l not in fermi_lines] 
    ylim = ax.get_xlim()
    ax.set_xlim()
    ax.set_ylim()
    ax.set_xlim(0,ymax)
    ax.set_ylim(ylim)
    ax.set_xlabel('Density of states',size=25)
    ax.set_yticks([])
    ax.set_ylabel(None)
    plt.xticks(fontsize=25)

    fig = plt.gcf()
    fig.tight_layout()

    return plt

