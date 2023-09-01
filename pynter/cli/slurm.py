#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:45:20 2023

@author: villa
"""

from pynter.slurm.core import Slurm
from pynter.slurm.job_settings import JobSettings

def setup_job_script(subparsers):
    
    parser_sub = subparsers.add_parser('job-script',help='Create job script')
    parser_sub = parse_job_script_args(parser_sub)
    parser_sub.set_defaults(func=write_job_script)
    return
    

def parse_job_script_args(parser):
    shdef = JobSettings()
    parser.add_argument('--filename',help='File name (default: %(default)s)',required=False,default=shdef.filename,type=str,metavar='',dest='filename')
    parser.add_argument('--array-size',help='Size of job array',required=False,default=shdef.array_size,type=int,metavar='',dest='array_size')
    parser.add_argument('--modules',action='append',help="Modules to load (default: %(default)s)" ,required=False,default=shdef.modules,type=str,metavar='',dest='modules')
    parser.add_argument('--exe',help='Path to executable (default: %(default)s)',required=False,default=shdef.path_exe,type=str,metavar='',dest='path_exe')
    parser.add_argument('--stop-array',action='store_true',help='Add lines to stop array jobs when converged (default: %(default)s)',required=False,default=False,dest='add_stop_array')
    parser.add_argument('--automation',help='Script with automation to add',required=False,default=shdef.add_automation,type=str,metavar='',dest='add_automation')
    parser.add_argument('--header',action='append',help='Add line to header part of script',required=False,default=shdef.add_lines_header,type=str,metavar='',dest='add_lines_header')
    parser.add_argument('--body',action='append',help='Add line to body part of script',required=False,default=shdef.add_lines_body,type=str,metavar='',dest='add_lines_body')

    slurm = Slurm()
    for key in slurm.arguments:
        if key in slurm.keys():
            parser.add_argument(f'--{key}',required=False,type=str,default=slurm[key],metavar='',dest=key)
        else:
            parser.add_argument(f'--{key}',required=False,type=str,metavar='',dest=key)

    return parser 


def write_job_script(args):
    kwargs = vars(args)
    if 'func' in kwargs.keys():
        del kwargs['func']

    JobSettings(**kwargs).write_bash_file()
    return
    