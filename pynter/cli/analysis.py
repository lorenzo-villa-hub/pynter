#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:07:47 2023

@author: villa
"""

import os
import warnings
import numpy as np

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.transition_state import NEBAnalysis

from pynter.tools.utils import save_object_as_json


def setup_analysis(subparsers):
    
    parser_analyse = subparsers.add_parser('analysis',help='Analysis tools for VASP')
    subparsers_analysis = parser_analyse.add_subparsers()
    
    parser_vasprun = subparsers_analysis.add_parser('vasprun',help='Analyse vasprun.xml file with pymatgen')
    setup_analyse_vasprun(parser_vasprun)
    
    parser_neb = subparsers_analysis.add_parser('NEB',help='Analyse NEB VASP calculation with pymatgen')
    setup_analyse_neb(parser_neb)
    

def setup_analyse_neb(parser_sub):    
    parser_sub.add_argument('-p','--path',help='Path to VASP calculation',required=False,type=str,default=None,metavar='',dest='path')
    parser_sub.set_defaults(func=analyse_neb)    
    return


def analyse_neb(args):
    path = args.path if args.path else os.getcwd()
    neb = NEBAnalysis.from_dir(path)
    print('Energies:\n',neb.energies)
    print('Forces:\n',neb.forces)
    return

    
def setup_analyse_vasprun(parser_sub):    
    parser_sub.add_argument('-p','--path',help='Path to VASP calculation',required=False,type=str,default=None,metavar='',dest='path')
    parser_sub.add_argument('--bandgap',help='Analyse band gap',required=False,default=False,action='store_true',dest='bandgap') 
    parser_sub.add_argument('--convergence',help='Analyse convergence of VASP calculation',required=False,default=False,action='store_true',dest='convergence') 
    parser_sub.add_argument('--dielectric-properties',help='Analyse dielectric properties',required=False,default=False,action='store_true',dest='dielectric_properties') 
    parser_sub.add_argument('--export-dos',help='Export DOS object as json file',required=False,default=False,action='store_true',dest='export_dos') 
    parser_sub.set_defaults(func=analyse_vasprun)
    return


def analyse_vasprun(args):
    path = args.path if args.path else os.getcwd()
    try:
        vasprun = Vasprun(os.path.join(path,'vasprun.xml'))
    except:
        vasprun = None
        warnings.warn('Vasprun could not be read, the calculation probably failed or is not finished\n',UserWarning)   
    if args.bandgap:
        analyse_bandgap(vasprun)
    if args.convergence:
        analyse_convergence(vasprun)
    if args.dielectric_properties:
        analyse_dielectric_properties(vasprun)
    if args.export_dos:
        dos = vasprun.complete_dos
        save_object_as_json(dos, 'complete_dos.json')
        
    return 
    
def analyse_bandgap(vasprun):
    gap, cbm, vbm, is_direct = vasprun.eigenvalue_band_properties
    gap = np.around(gap,decimals=4)
    print(f'Energy gap is {gap} eV')
    if is_direct:
      print('Direct gap')
    else:
      print('Indirect gap')
    print(f'CBM is at {cbm} eV')
    print(f'VBM is at {vbm} eV')
    return

def analyse_convergence(vasprun):
    if vasprun is None:
        conv_el = False
        conv_ionic = False
    else:
    	conv_el = vasprun.converged_electronic
    	conv_ionic = vasprun.converged_ionic
    print(f'Electronic convergence: {conv_el}')
    print(f'Ionic convergence: {conv_ionic}')
    return
    
def analyse_dielectric_properties(vasprun):
    incar = vasprun.incar
    print('')
    if 'LOPTICS' in incar:
        if incar['LOPTICS'] == True:
            print('"LOPTICS = .TRUE." found in INCAR \n')
            epsilon = vasprun.dielectric
            real_part_0 = epsilon[1][0]
            print('Real part of dielectric function at 0 frequency:')
            print(real_part_0)
            print('')
    if 'LEPSILON' in incar:
        if incar['LEPSILON'] == True:        
            print('"LEPSILON = .TRUE." found in INCAR \n')
            if incar['IBRION'] == 8:
                print('"IBRION = 8" found in INCAR \n')
            print('Macroscopic static dielectric tensor:')
            print(vasprun.epsilon_static)
            print('ionic part of the static dielectric constant. Present when it’s a DFPT run:')
            print(vasprun.epsilon_ionic)
            print('')
    return
