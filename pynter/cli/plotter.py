#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 14:04:04 2023

@author: villa
"""

import os
import json
import argparse

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter, BSDOSPlotter # real pymatgen 
from pymatgen.electronic_structure.dos import CompleteDos


#parser_sub = argparse.Argumentparser_sub(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

def setup_plot_dos(subparsers):

    parser_sub = subparsers.add_parser('plot-dos',help='Plot density of states with pymatgen')
    
    parser_sub.add_argument('-sg','--sigma',help='Standard deviation for Gaussian smearing of DOS',required=False,type=float,default=None,metavar='',dest='sigma')    
    parser_sub.add_argument('-s','--save',help='Save fig as pdf',required=False,default=False,action='store_true',dest='savefig')   
    parser_sub.add_argument('-st','--stack',help='Stack DOS',required=False,default=False,action='store_true',dest='stack') 
    parser_sub.add_argument('-y','--ylim',help='Range for y-axis',required=False,default=None,nargs='+',type=float,metavar='',dest='ylim')
    parser_sub.add_argument('-x','--xlim',help='Range for x-axis',required=False,default=None,nargs='+',type=float,metavar='',dest='xlim')
    parser_sub.add_argument('-p','--projection',help='Projection of DOS, "elements", "orbitals" or "None"',required=False,type=str,default='elements',metavar='',dest='dos_projection')    
    parser_sub.add_argument('-t','--title',help='Title of the plot',required=False,type=str,default=None,metavar='',dest='title') 
    parser_sub.set_defaults(func=plot_dos)
    
    return


def plot_dos(args):

    xlim = args.xlim if args.xlim else (-10,10)
    complete_dos = Vasprun('vasprun.xml').complete_dos
    
    plotter = DosPlotter(stack=args.stack,sigma=args.sigma)
    plotter.add_dos('Total',complete_dos)
    
    if args.dos_projection == 'orbitals':
        dos_dict = complete_dos.get_spd_dos()
    elif args.dos_projection == 'elements':
        dos_dict = complete_dos.get_element_dos()
    for k in dos_dict:
        if not args.dos_projection=='None':   
            plotter.add_dos(k,dos_dict[k])
    plt = plotter.get_plot(xlim=xlim,ylim=args.ylim)
    plt.title(args.title)
    if args.savefig:
        plt.savefig('DOS.pdf',bbox_inches='tight')
    else:
        plt.show()
        
        
def setup_plot_dos_bs(subparsers):

    parser_sub = subparsers.add_parser('plot-dos-bs',help='Plot density of states and band structure with pymatgen')
    
    parser_sub.add_argument('-p','--path-to-dos',help='Path to Dos object saved as json or path to where vasprun for dos is stored',required=False,type=str,default=None,metavar='',dest='path_to_dos')    
    parser_sub.add_argument('-s','--save',help='Save fig as pdf',required=False,default=False,action='store_true',dest='savefig')   
    parser_sub.add_argument('-hy','--hybrid-mode',help='Force hybrid mode for BS',required=False,default=False,action='store_true',dest='hybrid')   
    parser_sub.add_argument('-y','--ylim',help='Range for y-axis',required=False,default=None,nargs='+',type=float,metavar='',dest='ylim')
    parser_sub.add_argument('-t','--title',help='Title of the plot',required=False,type=str,default=None,metavar='',dest='title') 
    parser_sub.set_defaults(func=plot_dos_bs)
    return


def plot_dos_bs(args):

    if not args.ylim:
        args.ylim = (-10,10)    
    get_dos = False
    if args.path_to_dos:
        if '.json' in args.path_to_dos:
            with open(args.path_to_dos) as f:
                dos = CompleteDos.from_dict(json.load(f))
        elif args.path_to_dos:
            try:
                dos = Vasprun(os.path.join(args.path_to_dos,"vasprun.xml")).complete_dos
            except:
                raise ValueError(f'Reading of dos from vasprun.xml in {args.path_to_dos} failed')
    else:
        get_dos = True
    
    v = Vasprun('vasprun.xml')    
    if get_dos:
        dos = v.complete_dos
    efermi = dos.efermi
    gap = dos.get_gap()
    ylim = args.ylim
    vb_range = -1*ylim[0] if ylim else 4  # adjust for nonsense args in BSDOSPlotter
    cb_range = ylim[1] - gap if ylim else 4  # adjust for nonsense args in BSDOSPlotter
    bs = v.get_band_structure(line_mode=True,force_hybrid_mode=args.hybrid,efermi=efermi)
    
    plt = BSDOSPlotter(bs_projection=None, dos_projection='elements', bs_legend=None, 
                        vb_energy_range=vb_range, cb_energy_range=cb_range).get_plot(bs,dos)
    ax = plt.gcf().get_axes()[0]
    ax.set_ylabel('$E - E_F$ (eV)')
    
    if args.savefig:
        plt.savefig('DOS-BS.pdf')
    else:
        plt.show()
    
    
    
    
    
    
    