#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 11:27:24 2023

@author: villa
"""

from pymatgen.io.vasp.inputs import Poscar


def setup_mp_database(subparsers):
    
    parser_sub = subparsers.add_parser('MP',help='Interface with Materials Project database')
    
    parser_sub.add_argument('-id','--mp-id',help='Materials project ID',required=True,type=str,metavar='',dest='mp_id')
    parser_sub.add_argument('-s','--structure',help='print Structure',action='store_true',required=False,default=False,dest='print_structure')
    parser_sub.add_argument('-P','--poscar',help='create POSCAR',action='store_true',required=False,default=False,dest='write_poscar')
    parser_sub.add_argument('-c','--conv-cell',help='Get conventional unit cell (default: %(default)s)',action='store_true',
                        required=False,default=False,dest='get_conventional')
    parser_sub.add_argument('-f','--final-structure',help='Get final structure (default: %(default)s)',action='store_true',
                        required=False,default=False,dest='get_final')     
    parser_sub.add_argument('-n','--filename',help='Filename for POSCAR (default: %(default)s)',required=False,
                        default='POSCAR',type=str,metavar='',dest='filename')
    
    parser_sub.set_defaults(func=run_mp_database)
    return
    

def run_mp_database(args):
    from pynter.tools.materials_project import MPDatabase
    mp = MPDatabase(mp_id=args.mp_id)
    structure = mp.get_structure(final=args.get_final,conventional_unit_cell=args.get_conventional)
    
    if args.print_structure:
        print(structure)
    if args.write_poscar:
        Poscar(structure).write_file(args.filename)
        
        