#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:45:20 2023

@author: villa
"""

from pynter.slurm.job_script import ScriptHandler

def setup_job_script(subparsers):
    
    parser_sub = subparsers.add_parser('job-script',help='Create job script')
    parser_sub = parse_job_script_args(parser_sub)
    parser_sub.set_defaults(func=write_job_script)
    return
    

def parse_job_script_args(parser):
    shdef = ScriptHandler() #default values
    parser.add_argument('-A','--project',help='Project ID (default: %(default)s)',required=False,default=shdef.project_id,type=str,metavar='',dest='project_id')
    parser.add_argument('-n','--name',help='Job name',required=False,default=shdef.name,type=str,metavar='',dest='name')
    parser.add_argument('-a','--array',help='Size of job array',required=False,default=shdef.array_size,type=int,metavar='',dest='array_size')
    parser.add_argument('-e','--email',help='Email address for job notification',required=False,default=shdef.email,type=str,metavar='',dest='email')
    parser.add_argument('-N','--nodes',help='Number of nodes (default: %(default)s)',required=False,default=shdef.nodes,type=int,metavar='',dest='nodes')
    parser.add_argument('-c','--cores-per-node',help='Number of cores per node (default: %(default)s)',required=False,default=shdef.cores_per_node,type=int,metavar='',dest='cores_per_node')
    parser.add_argument('-out','--output',help='Output filename',required=False,default=shdef.output_filename,type=str,metavar='',dest='output_filename')
    parser.add_argument('-err','--error',help='Error filename',required=False,default=shdef.error_filename,type=str,metavar='',dest='error_filename')
    parser.add_argument('-t','--timelimit',help='Timelimit (default: %(default)s)',required=False,default=shdef.timelimit,type=str,metavar='',dest='timelimit')
    parser.add_argument('-M','--memory-per-cpu',help='Memory per cpu (default: %(default)s)',required=False,default=shdef.memory_per_cpu,type=int,metavar='',dest='memory_per_cpu')
    parser.add_argument('-p','--partition',help='Partition (default: %(default)s)',required=False,default=shdef.partition,type=str,metavar='',dest='partition')
    parser.add_argument('-C','--processor',help='avx or avx2 (default: %(default)s)',required=False,default=shdef.processor,type=str,metavar='',dest='processor')
    parser.add_argument('-ml','--modules',action='append',help="Modules to load (default: %(default)s)" ,required=False,default=shdef.modules,type=str,metavar='',dest='modules')
    parser.add_argument('-x','--exe',help='Path to executable (default: %(default)s)',required=False,default=shdef.path_exe,type=str,metavar='',dest='path_exe')
    parser.add_argument('-s','--stop-array',action='store_true',help='Add lines to stop array jobs when converged (default: %(default)s)',required=False,default=False,dest='add_stop_array')
    parser.add_argument('-S','--automation',help='Script with automation to add',required=False,default=shdef.add_automation,type=str,metavar='',dest='add_automation')
    parser.add_argument('-H','--header',action='append',help='Add line to header part of script',required=False,default=shdef.add_lines_header,type=str,metavar='',dest='add_lines_header')
    parser.add_argument('-B','--body',action='append',help='Add line to body part of script',required=False,default=shdef.add_lines_body,type=str,metavar='',dest='add_lines_body')
    parser.add_argument('-f','--filename',help='File name (default: %(default)s)',required=False,default=shdef.filename,type=str,metavar='',dest='filename')
    return parser


def write_job_script(args):
    kwargs = vars(args)
    if 'func' in kwargs.keys():
        del kwargs['func']
    sh = ScriptHandler(**kwargs)
    sh.write_script()
    return
    