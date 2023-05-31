#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 14:04:04 2023

@author: villa
"""

import os
import json

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter, BSDOSPlotter, BSPlotter 
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.analysis.transition_state import NEBAnalysis


def setup_plotter(subparsers):
    
    parser_plot = subparsers.add_parser('plot',help='Plotting tools for DOS, BS and NEB')
    subparsers_plot = parser_plot.add_subparsers()
    
    parser_bs = subparsers_plot.add_parser('BS',help='Plot band structure with pymatgen')
    setup_plot_bs(parser_bs)
    
    parser_dos = subparsers_plot.add_parser('DOS',help='Plot density of states with pymatgen')
    setup_plot_dos(parser_dos)
    
    parser_dos_bs = subparsers_plot.add_parser('DOS-BS',help='Plot density of states and band structure with pymatgen')
    setup_plot_dos_bs(parser_dos_bs)
    
    parser_neb = subparsers_plot.add_parser('NEB',help='Plot NEB energy landscape with pymatgen')
    setup_plot_neb(parser_neb)
    
    return
    



def setup_plot_bs(parser_sub):    
    parser_sub.add_argument('-k','--kpoints-filename',help='KPOINTS filename',required=False,type=str,default=None,metavar='',dest='kpoints_filename')    
    parser_sub.add_argument('-s','--save',help='Save fig as pdf',required=False,default=False,action='store_true',dest='savefig')   
    parser_sub.add_argument('-hy','--hybrid-mode',help='Force hybrid mode for BS',required=False,default=False,action='store_true',dest='force_hybrid_mode') 
    parser_sub.add_argument('-l','--line-mode',help='Line mode for BS',required=False,default=True,action='store_true',dest='line_mode')   
    parser_sub.add_argument('-y','--ylim',help='Range for y-axis',required=False,default=None,nargs='+',type=float,metavar='',dest='ylim')
    parser_sub.add_argument('-sm','--smooth',help='Smooth band structure',required=False,default=False,action='store_true',dest='smooth')  
    parser_sub.add_argument('-t','--title',help='Title of the plot',required=False,type=str,default=None,metavar='',dest='title') 
    parser_sub.set_defaults(func=plot_bs) 
    return


def plot_bs(args):
    v = Vasprun("vasprun.xml")
    bs = v.get_band_structure(kpoints_filename=args.kpoints_filename,line_mode=args.line_mode,
                              force_hybrid_mode=args.force_hybrid_mode)
    plt = BSPlotter(bs).get_plot(ylim=args.ylim,smooth=args.smooth)
    plt.gca().get_legend().remove()
    
    if args.title:
        plt.title(args.title,fontdict={'fontsize':25})
    if args.savefig:
        plt.savefig('BS.pdf',bbox_inches='tight')
    else:
        plt.show()
    return    
    




def setup_plot_dos(parser_sub): 
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
        
 
    
 
        
def setup_plot_dos_bs(parser_sub):    
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
    
    
    
    
    
def setup_plot_neb(parser_sub):   
    parser_sub.add_argument('-p','--path',help='Path to NEB calculation',required=False,type=str,default=None,metavar='',dest='path')
    parser_sub.add_argument('-pr','--print',help='Print energies and forces',required=False,default=False,action='store_true',dest='print_energies')
    parser_sub.add_argument('-s','--save',help='Save fig as pdf',required=False,default=False,action='store_true',dest='savefig')
    parser_sub.add_argument('-t','--title',help='Title of the plot',required=False,type=str,default=None,metavar='',dest='title')
    parser_sub.set_defaults(func=plot_neb)
    return
    
    
def plot_neb(args):
    path = args.path if args.path else os.getcwd()
    path = os.path.abspath(path)
    neb_analysis = NEBAnalysis.from_dir(path)
    if args.print_energies:
        print('Energies:\n',neb_analysis.energies)
        print('Forces:\n',neb_analysis.forces)
    plt = neb_analysis.get_plot()
    plt.title(args.title)
    if args.savefig:
        plt.savefig('NEB.pdf')
    else:
        plt.show()
        
    
    
    
    
    
    
    
    
    