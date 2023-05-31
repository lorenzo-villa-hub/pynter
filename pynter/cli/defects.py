#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:09:45 2023

@author: villa
"""

import os
from glob import glob
import numpy as np
import warnings

from pynter.data.datasets import Dataset
from pynter.vasp.jobs import VaspJob
from pynter.slurm.job_script import ScriptHandler
from pynter.defects.entries import DefectEntry
from pynter.defects.analysis import DefectsAnalysis
from pynter.defects.corrections import get_kumagai_correction_from_jobs

def setup_defects(subparsers):
    
    parser_defects = subparsers.add_parser('defects',help='Import and analyse defect calculations. Use with extreme care.')
    subparsers_defects = parser_defects.add_subparsers()
    
    parser_import = subparsers_defects.add_parser('import',help='Create defect entries from VASP DFT calculations')
    setup_import(parser_import)
    return
    




def setup_import(parser):
    job_script_filename = ScriptHandler().filename
    parser.add_argument('-pb','--path-bulk',help='Path to bulk calculation',required=True,type=str,metavar='',dest='path_bulk')
    parser.add_argument('-p','--path',help='Path to defect calculations, can contain wildcards (default: %(default)s)',required=False,type=str,default=os.getcwd(),metavar='',dest='path')
    
    parser.add_argument('-e','--exclude',action='append',help='Exclude specific defect types (Vacancy, Substitution, Interstitial, Polaron, DefectComplex)',
                        required=False,default=None,dest='exclude')
    parser.add_argument('-c','--corrections',action='store_true',help='Compute Kumagai corrections (default: %(default)s)',required=False,default=False,dest='corrections')
    parser.add_argument('-dt','--dielectric-tensor',help='Dielectric tensor, write the matrix in a line (a11 a12 a13 a21 a22 a23 a31 a32 a33)',
                        required=False,type=str,default=None,metavar='',dest='dielectric_tensor')    
    parser.add_argument('-j','--job-script-filename',help='Job script filename (default: %(default)s)',required=False,type=str,default=job_script_filename,metavar='',
                        dest='job_script_filename')
    parser.add_argument('-t','--tolerance',help='Tolerance in A° for automatic defect finding (default: %(default)s)',required=False,type=float,default=0.01,metavar='',
                        dest='tolerance')

    parser.set_defaults(func=import_entries)
    return


def import_entries(args):
    job_bulk = VaspJob.from_directory(args.path_bulk)
    path = args.path
    script = args.job_script_filename
    
    if '*' in path:
        dirs = glob(os.path.normpath(path))
        jobs = [VaspJob.from_directory(path=jdir,job_script_filename=script) for jdir in dirs]
        ds = Dataset(jobs)
    else:
        ds = Dataset.from_directory(path,script)
    print(ds)
    def get_last_node(job):
        if len(job.node_points) == 0:
            return None
        else:
            return job.node_points[-1]
    jobs_to_import = []
    for group in ds.groups:
        jobs = ds.select_jobs(groups=[group],is_converged=True)
        charges = []
        for j in jobs:
            if j.charge not in charges:
                charges.append(j.charge)
        for q in charges:
            jobs_same_charge = ds.select_jobs(jobs,charge=q)
            sorted_jobs_same_charge = ds.sort_jobs(jobs_same_charge,features=get_last_node)
            job_defect = sorted_jobs_same_charge[-1]
            if job_defect.formula != job_bulk.formula or job_defect.charge != job_bulk.charge:
                jobs_to_import.append(job_defect)
    
    
    entries = []
    tol = args.tolerance
    for job_defect in jobs_to_import:
        print(job_defect)
        corrections = {}
        if args.corrections:
            if args.dielectric_tensor:
                dielectric_tensor_string = ' '.join(args.dielectric_tensor.split()) # remove extra spaces
                dielectric_tensor = np.fromstring(dielectric_tensor_string, dtype=float, sep=' ').reshape((3, 3))
            else:
                raise ValueError('Dielectric tensor needs to be provided to perform Kumagai corrections')
            kumagai_corrections = get_kumagai_correction_from_jobs(job_defect, job_bulk, dielectric_tensor,tol=tol)
            kumagai_total = sum([v for v in kumagai_corrections.values()])
            corrections['kumagai'] = kumagai_total
        data = {'stress':job_defect.stress}
        entry = DefectEntry.from_jobs(job_defect,job_bulk,corrections=corrections,data=data,multiplicity=None,tol=tol)
        
        if entry.defect.defect_type == 'DefectComplex' and len(entry.defect.defects) > 3:
            warnings.warn('A complex with more than 3 species has been found, something has likely gone wrong in the automatic defect determination. Excluding it for now...')
        elif args.exclude==None or entry.defect.defect_type not in args.exclude:
            entries.append(entry)
    
    band_gap, vbm = job_bulk.energy_gap, job_bulk.vbm
    da = DefectsAnalysis(entries, vbm, band_gap)
    print(da)
    return
            
        

    
    