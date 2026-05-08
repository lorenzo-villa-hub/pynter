#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 12:38:27 2021

@author: lorenzo
"""
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt


def plot_NEB(
        neb_analysis,
        ax=None,
        normalize_rxn_coordinate: bool = True,
        label_barrier: bool = True,
        linewidth=2,
        markersize=10,
        fontsize=None,
        **kwargs):
    """
    Get NEB plot. Adapted code from ``pymatgen.analysis.transition_state.NEBAnalysis.get_plot``.

    Parameters
    ----------
    neb_analysis : NEBAnalysis
        NEBAnalysis object
    ax : Axes
        Matplotlib Axes
    normalize_rxn_coordinate : bool 
        Whether to normalize the reaction coordinate to between 0 and 1. Defaults to True.
    label_barrier : bool
        Whether to label the maximum barrier. Defaults to True.
    linewith : float
        Line width for spline.
    Markersize : float
        Marker size for energy points
    fontsize : float
        Font size.

    Returns
    -------
    matplotlib Axes object.
    """

    ax= _get_ax(ax=ax,**kwargs)
    _style_ax(ax,fontsize=fontsize,legend=False,grid=False)

    neb = neb_analysis
    scale = 1 / neb.r[-1] if normalize_rxn_coordinate else 1
    xs = np.arange(0, np.max(neb.r), 0.01)
    ys = neb.spline(xs)
    relative_energies = neb.energies - neb.energies[0]
    ax.plot(neb.r * scale, relative_energies, "ro", xs * scale, ys, "k-", linewidth=linewidth, markersize=markersize)

    ax.set_xlabel("Reaction Coordinate")
    ax.set_ylabel("Energy (meV)")
    ax.set_ylim((np.min(ys) - 0.01, np.max(ys) * 1.02 + 0.02))

    ax.set_title(f"$\\Delta E = $ {np.max(ys) - np.min(ys):.3f} eV")

    return ax  




def update_fontsize(size):
    matplotlib.rcParams.update({'font.size': size})
    return


def _get_ax(ax=None, **fig_kwargs):
    if ax is None:
        _, ax = plt.subplots(**fig_kwargs)
    return ax

def _style_ax(ax, fontsize=None, legend=True, grid=True):

    if fontsize is not None:
        ax.set_xlabel(ax.get_xlabel(), fontsize=fontsize)
        ax.set_ylabel(ax.get_ylabel(), fontsize=fontsize)
        ax.tick_params(axis='both', labelsize=fontsize*0.9)
        if legend:
            ax.legend(fontsize=fontsize*0.9)
    else:
        if legend:
            ax.legend()
    if grid:
        ax.grid()
    return