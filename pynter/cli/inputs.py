#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:05:08 2023

@author: villa
"""
import os
import importlib

from pynter.slurm.job_script import ScriptHandler
from pynter.tools.materials_project import MPDatabase
from pynter.vasp.schemes import InputSets,Schemes


def setup_inputs(subparsers):
    
    parser_schemes = subparsers.add_parser('inputs',help='Create inputs for several calculation schemes')
    subparsers_schemes = parser_schemes.add_subparsers()

    parser_vasp = subparsers_schemes.add_parser('vasp',help='Create inputs for VASP DFT calculations')
    setup_inputs_vasp(parser_vasp)
    
    return


def setup_inputs_vasp(parser):
    parser = parse_common_args(parser)
    parser = parse_vasp_args(parser)
    parser.add_argument('-is','--inputsets',help='Name of a set of inputs from the InputSets class',required=False,type=str,
                        default=None,metavar='',dest='input_sets')
    parser.set_defaults(func=create_inputs_vasp)
    return

def create_inputs_vasp(args):
    if args.poscar:
        structure = args.poscar.structure
    elif args.mp_id:
        structure = MPDatabase(args.mp_id).get_structure()
    else:
        raise ValueError('Either POSCAR or Materials Project ID needs to be provided')
    if args.job_script:
        script_path = os.path.abspath(os.path.dirname(args.job_script))
        job_settings = ScriptHandler.from_file(path=script_path,filename=os.path.basename(args.job_script)).settings
    else:
        job_settings = None
    
    if args.input_sets:
        input_sets = InputSets(path=args.path,structure=structure,incar_settings=args.incar,
                               kpoints=args.kpoints,potcar=args.potcar,job_settings=job_settings,
                               name=args.name,add_parent_folder=True)
        method = getattr(input_sets, args.input_sets)
        job = method()
        job.write_input()
    else:
        raise ValueError('Either method of InputSets or scheme must be provided')
    return
        
    


def parse_common_args(parser):
    parser.add_argument('-id','--materials-project-id',help='Materials project ID containing the structure, in case the POSCAR is not provided (default %(default)s)',
                        required=False,type=str,default=None,metavar='',dest='mp_id')

    parser.add_argument('-n','--name',help='Name of the folder containing input files (default %(default)s)',required=False,type=str,
                        default='inputs',metavar='',dest='name')
    
    parser.add_argument('-p','--path',help='Path for input files (default %(default)s)',required=False,type=str,
                        default=os.getcwd(),metavar='',dest='path')
    
    return parser


def parse_vasp_args(parser):

    parser.add_argument('-in','--incar',help='INCAR file with general settings, if not provided default settings are used',required=False,type=str,default=None,metavar='',
                        dest='incar')
    
    parser.add_argument('-kp','--kpoints',help='KPOINTS file, if not provided default settings are used',required=False,type=str,default=None,metavar='',
                        dest='kpoints')
    
    parser.add_argument('-pos','--poscar',help='POSCAR file (default: %(default)s)',required=False,type=str,default=None,metavar='',
                        dest='poscar')
    
    parser.add_argument('-pot','--potcar',help='POTCAR file (default: %(default)s)',required=False,type=str,default=None,metavar='',
                        dest='potcar')

    parser.add_argument('-j','--job-script',help='Job script with input job settings',required=False,type=str,
                        default=None,metavar='',dest='job_script')
    
    return parser  
    


