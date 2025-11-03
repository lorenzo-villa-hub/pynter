#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:30:49 2023

@author: villa
"""
import os

from pynter import SETTINGS
from pynter.jobs.vasp.vasp_automations import VaspAutomation
from pynter.jobs.datasets import Dataset

def parse_common_args(parser):
    auto = VaspAutomation()
    
    parser.add_argument('-j','--job-script',help='Job script filename (default: %(default)s)',required=False,default=SETTINGS['job_script_filename'],type=str,metavar='',dest='job_script_filename')
    parser.add_argument('-s','--status',help='Write exit status to file (default: %(default)s)',required=False,default=auto.status_filename,type=str,metavar='',dest='status_filename')
    parser.add_argument('-e','--error-check',action='store_true',help='Perform error checking (default: %(default)s)',required=False,default=False,dest='error_check')
    return parser


def setup_automation(subparsers):
    
    subparsers_automation = subparsers.add_parser('automation',help='Set up automation for calculations')
    subparsers_automation = subparsers_automation.add_subparsers()
    
    parser_vasp = subparsers_automation.add_parser('vasp',help='VASP automation')
    parser_vasp = parse_common_args(parser_vasp)
    parser_vasp.add_argument('-c','--contcar',action='store_true',help='Copy CONTCAR to POSCAR of next step',required=False,default=False,dest='contcar')
    parser_vasp.add_argument('-W','--wavecar',action='store_true',help='Copy WAVECAR to next step',required=False,default=False,dest='wavecar')
    parser_vasp.add_argument('-C','--chgcar',action='store_true',help='Copy CHGCAR to next step',required=False,default=False,dest='chgcar')
    parser_vasp.add_argument('-K','--check-kpoints',action='store_true',help='Copy WAVECAR and POSCAR only if KPOINTS of next step are the same',required=False,default=False,dest='check_kpoints')
    parser_vasp.set_defaults(func=run_automation_vasp)
    
    parser_vasp_neb = subparsers_automation.add_parser('vasp-NEB',help='VASP automation for NEB calculations')
    parser_vasp_neb = parse_common_args(parser_vasp_neb)
    parser_vasp_neb.add_argument('-c','--contcar',action='store_true',help='Copy CONTCAR to POSCAR of next step',required=False,default=False,dest='contcar')
    parser_vasp_neb.add_argument('-W','--wavecar',action='store_true',help='Copy WAVECAR to next step',required=False,default=False,dest='wavecar')
    parser_vasp_neb.add_argument('-C','--chgcar',action='store_true',help='Copy CHGCAR to next step',required=False,default=False,dest='chgcar')
    parser_vasp_neb.add_argument('-K','--check-kpoints',action='store_true',help='Copy WAVECAR and POSCAR only if KPOINTS of next step are the same',required=False,default=False,dest='check_kpoints')
    parser_vasp_neb.set_defaults(func=run_automation_vasp_neb)


def run_automation_vasp(args):
    automation = VaspAutomation()
    ds = Dataset.from_directory(
                                path='../',
                                job_script_filenames=args.job_script_filename,
                                sort='path')
    path = os.getcwd()
    job_current = ds.select_jobs(path=path)[0]
    is_converged_electronic, is_converged_ionic = automation.check_convergence(job_current)
    if is_converged_electronic and is_converged_ionic:
        current_index = ds.jobs.index(job_current)
        job_next = ds.jobs[current_index + 1]
        files = []
        if args.contcar:
            files.append('CONTCAR')
        if args.chgcar:
            files.append('CHGCAR')
        if args.wavecar:
            files.append('WAVECAR')
        automation.copy_files_from_job1_to_job2(
                                    job1=job_current,
                                    job2=job_next,
                                    filenames=files,
                                    check_kpoints=args.check_kpoints)
        
        job_next.run_job(sync=False,write_input=False)
    automation.write_status(path=path)


def run_automation_vasp_neb(args):
    raise NotImplementedError('VASP NEB automation still to be implemented')