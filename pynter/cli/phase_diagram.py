#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 14:50:11 2023

@author: villa
"""
import os

from pymatgen.core.composition import Composition
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram

from pynter import SETTINGS
from pynter.phase_diagram.chempots import PDHandler, Reservoirs
from pynter.tools.utils import save_object_as_json


def setup_phase_diagram(subparsers):
    
    parser_sub = subparsers.add_parser('phase-diagram',help='Generate phase diagrams with pymatgen from the Materials Project database')
    
    parser_sub.add_argument('-e','--elements',action='append',help='Components of the phase diagram, provide as "el1-el2" or -e el1 -e el2',
                            type=str,metavar='',dest='elements')
    parser_sub.add_argument('-c','--chempots',help='Print boundary chempots of target phase',required=False,default=False,type=str,metavar='',
                            dest='composition_chempots')
    parser_sub.add_argument('-sc','--save-chempots',action='store_true',help='Save chempots in Reservoirs object as json',required=False,default=False,dest='savechempots')
    parser_sub.add_argument('-p','--plot-PD',action='store_true',help='Plot phase diagram',required=False,default=False,dest='plot_pd')
    parser_sub.add_argument('-ps','--plot-stability',action='append',help='Plot stability diagram, provide elements to be on axis as "el1-el2" or -ps el1 -ps el2'
                            ,required=False,type=str,dest='elements_stability')
    parser_sub.add_argument('-j','--json',action='store_true',help='Save PhaseDiagram as json',required=False,default=False,dest='savejson')
    parser_sub.add_argument('-s','--savefig',action='store_true',help='Save matplotlib figure as pdf',required=False,default=False,dest='savefig')
    parser_sub.set_defaults(func=run_phase_diagram)
    return


def run_phase_diagram(args):
    elements = args.elements[0].split('-') if len(args.elements)==1 else args.elements
    plt = None
    pd = get_phase_diagram(elements)
    print('\n',pd)
    
    if args.composition_chempots:
        composition = Composition(args.composition_chempots)
        chempots = get_chempots(pd,composition)
        reservoirs = Reservoirs(chempots,pd,are_chempots_delta=False)
        print('\n',f'Chemical potentials of {args.composition_chempots}:\n',reservoirs)
        if args.savechempots:
            fname = 'reservoirs_boundary_' + args.composition_chempots + '.json'
            path = os.path.join(os.getcwd(),fname)
            reservoirs.to_json(path)
    
    if args.plot_pd:
        plt = plot_pd(pd)
        if args.savefig:
            fname = 'PD_' + '-'.join(elements) + '.pdf'
            path = os.path.join(os.getcwd(),fname)
            plt.savefig(path,bbox_inches='tight')
        else:
            plt.show()
    if args.elements_stability:
        elements_stability = args.elements_stability[0].split('-') if len(args.elements_stability)==1 else args.elements_stability
        plt_stab = plot_stability_diagram(pd,elements_stability)
        if args.savefig:
            fname = 'SD_' + '-'.join(elements) + '.pdf'
            path = os.path.join(os.getcwd(),fname)
            plt_stab.savefig(path,bbox_inches='tight')
        else:
            plt_stab.show()
    
    if args.savejson:
        fname = 'PD_' + '-'.join(elements) + '.json'
        path = os.path.join(os.getcwd(),fname)
        save_object_as_json(pd,path)
    
    return


def get_chempots(pd,composition):
    return PDHandler(pd).get_all_boundaries_chempots(composition)
    
def get_phase_diagram(elements):  
    from pymatgen.ext.matproj import MPRester
    with MPRester(SETTINGS['API_KEY']) as mpr:
        compat = MaterialsProjectCompatibility()  # sets energy corrections and +U/pseudopotential choice
        unprocessed_entries = mpr.get_entries_in_chemsys(elements,inc_structure=True)
    processed_entries = compat.process_entries(unprocessed_entries)  # filter and add energy corrections    
    return PhaseDiagram(processed_entries)

def plot_pd(pd):
    return PDHandler(pd).get_plot()

def plot_stability_diagram(pd,axis_elements):
    return PDHandler(pd).get_stability_diagram(axis_elements)


   