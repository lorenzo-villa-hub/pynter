#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:07:47 2023

@author: villa
"""

import os
import warnings
import numpy as np
from pathlib import Path

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.transition_state import NEBAnalysis

from pynter.tools.utils import save_object_as_json


def setup_parse(subparsers):
    
    parser_parse = subparsers.add_parser('parse',help='Parse calculations')
    subparsers_parse = parser_parse.add_subparsers()
    
    parser_vasp = subparsers_parse.add_parser('vasp',help='Parse VASP calculations')
    setup_parser_vasp(parser_vasp)
    
    #parser_lammps = subparsers_parse.add_parser('LAMMPS',help='Parse LAMMPS calculations')
    #setup_parser_lammps(parser_lammps)
    

def setup_parser_vasp(parser_sub):    
    parser_sub.add_argument('-p','--path',help='Path to VASP calculations',required=False,type=str,default=None,metavar='',dest='path')
    parser_sub.add_argument('--ase',help='Store VASP calculations in ASE database',required=False,default=False,action='store_true',dest='ase_db') 
    parser_sub.set_defaults(func=parse_vasp_dirs)    
    return



def parse_vasp_dirs(args):
    if args.ase_db:
        import ase.db
        import ase.io
        target_path = Path(args.path) if args.path else Path(os.getcwd())
        db_name = f'{target_path.name}.db'
        print(f'Connecting to database: {db_name}')
        db = ase.db.connect(db_name)
        for path in target_path.rglob('*vasprun.xml'):
            print(f'Reading calculation in {path}')
            atoms = ase.io.read(str(path),index=-1)
            db.write(atoms=atoms,path=str(path))
    return


    
# def setup_analyse_vasprun(parser_sub):    
#     parser_sub.add_argument('-p','--path',help='Path to VASP calculation',required=False,type=str,default=None,metavar='',dest='path')
#     parser_sub.add_argument('--bandgap',help='Analyse band gap',required=False,default=False,action='store_true',dest='bandgap') 
#     parser_sub.add_argument('--convergence',help='Analyse convergence of VASP calculation',required=False,default=False,action='store_true',dest='convergence') 
#     parser_sub.add_argument('--dielectric-properties',help='Analyse dielectric properties',required=False,default=False,action='store_true',dest='dielectric_properties') 
#     parser_sub.add_argument('--export-dos',help='Export DOS object as json file',required=False,default=False,action='store_true',dest='export_dos') 
#     parser_sub.set_defaults(func=analyse_vasprun)
#     return
