#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:45:20 2023

@author: villa
"""

from pynter.hpc.slurm import Sbatch, JobSettings


def setup_job_script(subparsers):
    
    parser_sub = subparsers.add_parser('job-script',help='Create job script')
    parser_sub = parse_job_script_args(parser_sub)
    parser_sub.set_defaults(func=write_job_script)
    return
    

def parse_job_script_args(parser):
    shdef = JobSettings()
    parser.add_argument('--filename',help='File name (default: %(default)s)',required=False,default=shdef.filename,type=str,metavar='',dest='filename')
    parser.add_argument('--lines',action='append',help='Add line to script',required=False,default=shdef.script_lines,type=str,metavar='',dest='script_lines')

    sbatch = Sbatch()
    for key in sbatch.arguments:
        if key in sbatch.keys():
            parser.add_argument(f'--{key}',required=False,type=str,default=sbatch[key],metavar='',dest=key)
        else:
            parser.add_argument(f'--{key}',required=False,type=str,metavar='',dest=key)

    return parser 


def write_job_script(args):
    kwargs = vars(args)
    if 'func' in kwargs.keys():
        del kwargs['func']

    JobSettings(**kwargs).write_bash_file()
    return
    