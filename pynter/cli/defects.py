#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:09:45 2023

@author: villa
"""

import os
import json
from glob import glob
import numpy as np
import warnings
from pprint import pprint

from pymatgen.electronic_structure.dos import FermiDos

from pynter.data.datasets import Dataset
from pynter.vasp.jobs import VaspJob
from pynter.slurm.job_script import ScriptHandler
from pynter.defects.entries import DefectEntry
from pynter.defects.analysis import DefectsAnalysis
from pynter.defects.corrections import get_kumagai_correction_from_jobs
from pynter.phase_diagram.chempots import Chempots, Reservoirs
from pynter.tools.utils import save_object_as_json, get_object_from_json
from pynter.cli import inputs
from pynter.cli.utils import get_dict_from_line_string, round_floats



def setup_defects(subparsers):
    
    parser_defects = subparsers.add_parser('defects',help='Import and analyse defect calculations. Use with extreme care.')
    subparsers_defects = parser_defects.add_subparsers()

    parser_inputs = subparsers_defects.add_parser('inputs',help='Create inputs for VASP DFT calculations')
    setup_inputs(parser_inputs)
    
    parser_import = subparsers_defects.add_parser('import',help='Create defect entries from VASP DFT calculations')
    setup_import(parser_import)
 
    parser_plot = subparsers_defects.add_parser('plot',help='Plot collection of defect entries')
    setup_plot(parser_plot)

    parser_analysis = subparsers_defects.add_parser('analysis',help='Analyse collection of defect entries')
    setup_analysis(parser_analysis)

    return
    


def setup_inputs(parser):
    subparsers = parser.add_subparsers()
    
    parser_vasp = subparsers.add_parser('vasp',help='VASP DFT calculations')
    
    parser_vasp = inputs.parse_common_args(parser_vasp)
    parser_vasp = inputs.parse_vasp_args(parser_vasp)

    parser_vasp.add_argument('-auto','--automation',action='store_true',help='Add default automation to job script',required=False,
                        default=False,dest='automation')    

    parser_vasp.add_argument('-ss','--supercell-size',action='append',help='Size of the supercell',type=int,required=False,
                        default=None,metavar='',dest='supercell_size')

    parser_vasp.add_argument('-rel','--relaxation-scheme',help='Relaxation scheme to use, choose between "default" (2-step PBE) and "gamma" (4-step PBE) (default: %(default)s)',
                        required=False,default='default',type=str,metavar='',dest='relaxation_scheme')
    
    parser_vasp.add_argument('-sub','--substitutions',help='Substituions inputs. Provide elements and charges as {"<new_el>-on-<old_el>":[q0,q1,q2]}',type=str,required=False,
                        default=None,metavar='',dest='substitutions')
    
    parser_vasp.add_argument('-vac','--vacancies',help='Vacancies inputs. Provide elements and charges as {"el":[q0,q1,q2]}',type=str,required=False,
                        default=None,metavar='',dest='vacancies')
    
    parser_vasp.set_defaults(func=create_vasp_inputs)
    return


def create_vasp_inputs(args):
    jobs = []
    schemes = inputs.get_schemes(args)
    if args.substitutions:
        elements_to_replace_with_charges = json.loads(args.substitutions)
        jobs_sub = schemes.substitutions_pbe_relaxation(elements_to_replace_with_charges=elements_to_replace_with_charges,
                                                    supercell_size=args.supercell_size,automation=args.automation,
                                                    locpot=True,rel_scheme=args.relaxation_scheme)
        jobs = jobs + jobs_sub
    if args.vacancies:
        elements_with_charges = json.loads(args.vacancies)
        jobs_vac = schemes.vacancies_pbe_relaxation(elements_with_charges=elements_with_charges,
                                                    supercell_size=args.supercell_size,automation=args.automation,
                                                    locpot=True,rel_scheme=args.relaxation_scheme)
        jobs = jobs + jobs_vac

    ds = Dataset(jobs)
    ds.write_jobs_input()
    return


def setup_import(parser):
    job_script_filename = ScriptHandler().filename
    parser.add_argument('-pb','--path-bulk',help='Path to bulk calculation',required=True,type=str,metavar='',dest='path_bulk')
    parser.add_argument('-p','--path',help='Path to defect calculations, can contain wildcards (default: %(default)s)',required=False,type=str,default=os.getcwd(),metavar='',dest='path')
    
    parser.add_argument('-e','--exclude',action='append',help='Exclude specific defect types (Vacancy, Substitution, Interstitial, Polaron, DefectComplex)',
                        required=False,default=None,metavar='',dest='exclude')
    parser.add_argument('-c','--corrections',action='store_true',help='Compute Kumagai corrections (default: %(default)s)',required=False,default=False,dest='corrections')
    parser.add_argument('-dt','--dielectric-tensor',help='Dielectric constant or tensor, if tensor write the matrix in a line (a11 a12 a13 a21 a22 a23 a31 a32 a33)',
                        required=False,type=str,default=None,metavar='',dest='dielectric_tensor')    
    parser.add_argument('-j','--job-script-filename',help='Job script filename (default: %(default)s)',required=False,type=str,default=job_script_filename,metavar='',
                        dest='job_script_filename')
    parser.add_argument('-t','--tolerance',help='Tolerance in A° for automatic defect finding (default: %(default)s)',required=False,type=float,default=0.01,metavar='',
                        dest='tolerance')
    parser.add_argument('-s','--savejson',action='store_true',help='Save DefectsAnalysis and bulk DOS objects as json (default: %(default)s)',required=False,default=False,dest='savejson')

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
    print('\nCalculations imported from:')    
    print(ds.jobs_table(jobs_to_import,display=['charge']))
    
    entries = []
    tol = args.tolerance
    for job_defect in jobs_to_import:
        corrections = {}
        if args.corrections:
            if args.dielectric_tensor:
                if len(args.dielectric_tensor.split()) > 1:
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
            entries.append(entry)
    
    band_gap, vbm = job_bulk.energy_gap, job_bulk.vbm
    da = DefectsAnalysis(entries, vbm, band_gap)
    da.filter_entries(inplace=True,exclude=True,types=args.exclude)
    print('\nDefectsAnalysis:')
    print(da.get_dataframe())
    if args.savejson:
        da.to_json('defects_analysis_%s.json' %job_bulk.formula)
        job_bulk.get_output_properties(data=['complete_dos'])
        dos = job_bulk.complete_dos
        save_object_as_json(dos,'./DOS_%s.json'%job_bulk.formula)
        print('\nDefectsAnalysis object saved as defects_analysis_%s.json' %job_bulk.formula)
        print('DOS object saved as DOS_%s.json' %job_bulk.formula)
    return
  



def parse_common_args(parser):
    parser.add_argument('-f','--file',help='Json file representation of a DefectAnalysis object',
                        required=True,type=str,default=None,metavar='',dest='file')
    
    parser.add_argument('-c','--chempots',help="""Chemical potentials, written in a line ( el1 -x el2 -y el3 -z ). 
                        If not provided, elemental chempots are taken from the Materials Project. Use with care.""",
                        required=False,type=str,default=None,metavar='',dest='chempots')
    
    parser.add_argument('-e','--exclude',action='append',help='Exclude specific defect types (Vacancy, Substitution, Interstitial, Polaron, DefectComplex)',
                        required=False,default=None,metavar='',dest='exclude')
    
    parser.add_argument('-ee','--exclude-elements',action='append',help='Exclude defects containing these elements',required=False,
                    default=None,metavar='',dest='exclude_elements')
    return parser    

def get_defects_analysis(args):
    if args.file:
        da = DefectsAnalysis.from_json(args.file)
        if args.exclude:
            da.filter_entries(inplace=True,exclude=True,types=args.exclude)
        if args.exclude_elements:
            da.filter_entries(inplace=True,exclude=True,elements=args.exclude_elements)
    else:
        raise ValueError('DefectsAnalysis object must be provided as json file')
    return da
        
def get_chempots(args):
    if args.chempots:
        chempots_dict = get_dict_from_line_string(args.chempots)
    else:
        from pynter.tools.materials_project import MPDatabase
        mpd = MPDatabase()
        compositions = args.da.elements
        entries_dict = mpd.get_entries_from_compositions(compositions,lowest_e_above_hull=True)
        chempots_dict = {el:entry.energy_per_atom for el,entry in entries_dict.items()} #get elemental chempots from MP
        warnings.warn('Chemical potentials have not been provided, elemental chemical potentials taken ' \
                      'from the Materials Project database. Use with care.')
        print('Chemical potentials:')
        pprint(round_floats(chempots_dict))
    return Chempots(chempots_dict)
    


def setup_plot(parser):
    parser = parse_common_args(parser)

    parser.add_argument('-b','--binding-energies',action='store_true',help='Plot binding energies',required=False,
                    default=False,dest='plot_binding_energies')

    parser.add_argument('-ctl','--charge-transition-levels',action='store_true',help='Plot formation energies, provide chempots',required=False,
                    default=False,dest='plot_charge_transition_levels')

    parser.add_argument('-eform','--formation-energies',action='store_true',help='Plot formation energies, provide chempots',required=False,
                    default=False,dest='plot_formation_energies')
    
    parser.add_argument('-r','--reservoirs',help='Save plots as pdf',type=str,required=False,
                    default=None,dest='path_reservoirs')

    parser.add_argument('-s','--savefig',action='store_true',help='Save plots as pdf',required=False,
                    default=False,dest='savefig')
    
    parser.set_defaults(func=plot)
    return

def plot(args):
    da = get_defects_analysis(args)
    args.da = da
    if args.plot_binding_energies:
        if len(da.select_entries(types=['DefectComplex'])) == 0:
            raise ValueError('No DefectComplex found. Check that the defect complexes are present in the entries')
        else:
            plt = da.plot_binding_energies()
            if args.savefig:
                plt.savefig('binding_energies.pdf',bbox_inches='tight')
            else:
                plt.show()
    
    if args.plot_charge_transition_levels:
        plt = da.plot_ctl()
        if args.savefig:
            plt.savefig('charge_transition_levels.pdf',bbox_inches='tight')
        else:
            plt.show()
    
    if args.plot_formation_energies:
        if args.path_reservoirs:
            if args.chempots:
                warnings.warn('Reservoirs provided, ignoring input chemical potentials...')
            reservoirs = Reservoirs.from_json(args.path_reservoirs)
            nres = len(reservoirs)
            if nres <= 3:
                nrows,ncolumns = 1, nres
            else:
                ncolumns = 3
                nrows = nres//3 if nres%3==0 else nres//3 +1
            plotsize = (0.83*ncolumns,1.5*nrows)
            index=0
            for r,mu in reservoirs.items():
                index += 1
                plt = da.plot(chemical_potentials=mu,title=r,get_subplot=True,
                        plotsize=plotsize,subplot_settings=[nrows,ncolumns,index])
        else:
            chempots = get_chempots(args)
            plt = da.plot(chempots)
        if args.savefig:
            plt.savefig('formation_energies.pdf',bbox_inches='tight')
        else:
            plt.show()       
    return
    




def setup_analysis(parser):
    parser = parse_common_args(parser)        
    parser.add_argument('-d','--dos',help='Bulk DOS json file',type=str,required=False,default=None,metavar='',dest='dos')
    parser.add_argument('-ef','--fermi',help='Input Fermi level  (default: %(default)s)',required=False,metavar='',default=0,dest='fermi_level')
    parser.add_argument('-T','--temperature',help='Input temperature  (default: %(default)s)',required=False,default=300,metavar='',dest='temperature')
    
    parser.add_argument('-cc','--carrier-concentrations',action='store_true',help='Compute carrier concentrations, provide Fermi level, chempots and DOS json file',required=False,
                    default=False,dest='carrier_concentrations')
    
    parser.add_argument('-ctl','--charge-transition-levels',action='store_true',help='Compute charge transition levels',required=False,
                    default=False,dest='charge_transition_levels')
    
    parser.add_argument('-dc','--defect-concentrations',action='store_true',help='Compute defect concentrations, provide Fermi level and chempots',required=False,
                    default=False,dest='defect_concentrations')
    
    parser.add_argument('-eform','--formation-energies',action='store_true',help='Compute formation energies, provide chemical potentials',required=False,
                        default=False,dest='formation_energies')
    
    parser.add_argument('-sef','--solve-fermi-level',action='store_true',help='Solve charge neutrality to get Fermi level, provide chempots and DOS json file',required=False,
                        default=False,dest='solve_fermi_level')
    parser.set_defaults(func=analysis)
    return

def analysis(args):
    da = get_defects_analysis(args)
    args.da = da
    ef = args.fermi_level
    
    chempots = get_chempots(args)
    temperature = args.temperature
    
    if args.charge_transition_levels:
        ctl = da.charge_transition_levels()
        print('\nCharge transition levels:')
        pprint(round_floats(ctl))
    
    if args.formation_energies:
        formation_energies = da.formation_energies(chempots,ef)
        print('\nFormation energies:')
        pprint(round_floats(formation_energies))
    if args.defect_concentrations:
        if args.defect_concentrations:
            print('\nDefect concentrations:')
            print(da.defect_concentrations(chempots,temperature,ef))
        
    dos = get_object_from_json(FermiDos,args.dos) if args.dos else None
    if args.carrier_concentrations:
        if dos:
            print('\nCarrier concentrations:')
            print(round_floats(da.carrier_concentrations(dos,temperature,ef)))
        else:
            raise ValueError('DOS json file needs to be provided to compute carrier concentrations')
    if args.solve_fermi_level:
        if dos and chempots:
            print('\nEquilibrium Fermi level:')
            print(round(da.solve_fermi_level(chempots,dos,temperature),4))
        else:
            raise ValueError('DOS json file needs to be provided to compute carrier concentrations')
            
    return
        
    
    
        

        
    



















    
    