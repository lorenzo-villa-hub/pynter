#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 11:19:30 2021

@author: villa
"""
from pymatgen.electronic_structure.plotter import BSPlotter,BSDOSPlotter,DosPlotter, BSPlotterProjected
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


def plot_bs_elt_projected(bs,**kwargs):
    """
    Get element projected BS color plot with pymatgen. Note that vasprun.xml needs to
    be parsed with parse_projected_eigen set to True.

    Parameters
    ----------
    bs : 
        BandStructureSymmLine object, most likely generaten from Vasprun or BSVasprun.
    **kwargs : (dict)
        Arguments for the get_elt_projected_plots_color function in BSPlotter in pymatgen.

    Returns
    -------
    plt : 
        Matplotlib object.
    """
    plotter = BSPlotterProjected(bs)
    plt = plotter.get_elt_projected_plots_color(**kwargs)
    
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
    plt = plot_dos_dict(complete_dos, d, **kwargs)
    
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
    plt = plot_dos_dict(complete_dos, d, **kwargs)
    
    return plt


