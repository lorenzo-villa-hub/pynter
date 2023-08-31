#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 15:05:08 2023

@author: villa
"""
import os

from pymatgen.io.vasp.inputs import Incar,Kpoints,Poscar,Potcar

from pynter.data.datasets import Dataset
from pynter.slurm.job_settings import JobSettings
from pynter.vasp.schemes import AdvancedSchemes


def setup_inputs(subparsers):
    
    parser_schemes = subparsers.add_parser('inputs',help='Create inputs for several calculation schemes')
    subparsers_schemes = parser_schemes.add_subparsers()

    parser_vasp = subparsers_schemes.add_parser('vasp',help='Create inputs for VASP DFT calculations')
    setup_inputs_vasp(parser_vasp)
    
    return


def setup_inputs_vasp(parser):
    parser = parse_common_args(parser)
    parser = parse_vasp_args(parser)
    
    parser.add_argument('-cs','--charge-states',action='append',help='Charged states. Provide a list of integer charge states',required=False,
                        type=int,default=None,metavar='',dest='charge_states')
    
    parser.add_argument('-dp','--dielectric-properties',action='store_true',help='Electronic and ionic parts of the dielectric tensor',required=False,
                        default=False,dest='dielectric_properties')

    parser.add_argument('-ec','--encut-convergence',action='store_true',help='Energy cutoff convergence',required=False,
                        default=False,dest='encut_convergence')
    
    parser.add_argument('-fc','--fractional-charge',action='store_true',help='Fractional charge linearity',required=False,
                        default=False,dest='fractional_charge')
    
    parser.add_argument('-kc','--kpoints-convergence',action='store_true',help='Kpoints convergence',required=False,
                        default=False,dest='kpoints_convergence')
        
    parser.add_argument('-hse','--hse-relaxation',action='store_true',help='5 steps atomic and volume relaxation with HSE',required=False,
                        default=False,dest='hse_relaxation')
    
    parser.add_argument('-pbe','--pbe-relaxation',action='store_true',help='3 steps atomic and volume relaxation with PBE',required=False,
                        default=False,dest='pbe_relaxation')
    
    parser.add_argument('-pbe-es','--pbe-electronic-structure',action='store_true',help='3 steps electronic structure calculations with PBE',required=False,
                        default=False,dest='pbe_electronic_structure')
    
    parser.add_argument('-is','--inputsets',help='Name of a set of inputs from the InputSets class',required=False,type=str,
                        default=None,metavar='',dest='input_sets')
    
    parser.set_defaults(func=create_inputs_vasp)
    return


def create_inputs_vasp(args):
    schemes = get_schemes(args)
    jobs = None
    if args.encut_convergence:
        jobs = schemes.convergence_encut()
    
    elif args.kpoints_convergence:
        jobs = schemes.convergence_kpoints()
        
    elif args.charge_states:
        jobs = schemes.charge_states(args.charge_states,locpot=True)
        
    elif args.dielectric_properties:
        jobs = schemes.dielectric_properties_complete()
        
    elif args.fractional_charge:
        jobs = schemes.fractional_charge_linearity()
        
    elif args.pbe_relaxation:
        jobs = schemes.pbe_vol_relaxation()
        
    elif args.pbe_electronic_structure:
        jobs = schemes.pbe_electronic_structure()
    
    elif args.hse_relaxation:
        jobs = schemes.hse_vol_relaxation()
    
    if jobs:
        ds = Dataset(jobs)
        ds.write_jobs_input()        
 
    elif args.input_sets:
        method = getattr(schemes, args.input_sets)
        job = method()
        job.write_input()
    else:
        raise ValueError('Either method of InputSets or scheme must be provided')
    return
        

def get_schemes(args):    
    if args.poscar:
        poscar = Poscar.from_file(os.path.abspath(args.poscar))
        structure = poscar.structure
    elif args.mp_id:
        from pynter.tools.materials_project import MPDatabase
        structure = MPDatabase(args.mp_id).get_structure()
    else:
        raise ValueError('Either POSCAR or Materials Project ID needs to be provided')

    incar = Incar.from_file(os.path.abspath(args.incar)) if args.incar else None
    kpoints = Kpoints.from_file(os.path.abspath(args.kpoints)) if args.kpoints else None
    potcar = Potcar.from_file(os.path.abspath(args.potcar)) if args.potcar else None
    if args.job_script:
        script_path = os.path.abspath(os.path.dirname(args.job_script))
        job_settings = JobSettings.from_bash_file(path=script_path,filename=os.path.basename(args.job_script))
    else:
        job_settings = None
    
    schemes = AdvancedSchemes(path=args.path,structure=structure,incar_settings=incar,
                               kpoints=kpoints,potcar=potcar,job_settings=job_settings,
                               name=args.name,add_parent_folder=True)
    return schemes


def parse_common_args(parser):
    parser.add_argument('-id','--materials-project-id',help='Materials project ID containing the structure, in case the POSCAR is not provided (default: %(default)s)',
                        required=False,type=str,default=None,metavar='',dest='mp_id')

    parser.add_argument('-n','--name',help='Name of the folder containing input files (default: %(default)s)',required=False,type=str,
                        default='inputs',metavar='',dest='name')
    
    parser.add_argument('-p','--path',help='Path for input files (default: %(default)s)',required=False,type=str,
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
    


