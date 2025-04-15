#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 19:20:55 2023

@author: lorenzo
"""
import numpy as np
import warnings

from pymatgen.io.vasp.inputs import Potcar


def get_charge_from_computed_entry(entry):
    potcar_symbols_bulk = [d['titel'].split(' ')[1] for d in entry.parameters['potcar_spec']]
    potcar = Potcar(potcar_symbols_bulk)
    charge = 0
    if 'parameters' in entry.data.keys():
        parameters = entry.data['parameters']
    else:
        raise ValueError('"parameters" need to be present in ComputedEntry data')
    if 'NELECT' in parameters.keys():
        nelect = parameters['NELECT']
        val = {}
        for p in potcar:
            val[p.element] = p.nelectrons
        neutral = sum([ val[el.symbol]*coeff 
                       for el,coeff in entry.structure.composition.items()])
        charge = neutral - nelect
    if not isinstance(charge,int):
        charge = np.around(charge,decimals=1)
        
    return charge


def get_convergence_from_outcar(file='OUTCAR'):
    
    try:
        with open(file, "r") as f:
            lines = f.readlines()

        is_converged_electronic, is_converged_ionic, is_job_finished = False,False,False
        n_electronic_steps = 0
        for line in lines:    
            if "reached required accuracy" in line:
                n_electronic_steps += 1
            if "stopping structural energy minimisation" in line:
                is_converged_ionic = True
            if "General timing and accounting informations for this job" in line:
                is_job_finished = True
        
        parameters = get_parameters_from_outcar(tags=['IBRION','NELM'],file=file)
        ibrion = int(parameters['IBRION'])
        nelm = int(parameters['NELM'])
        
        if n_electronic_steps < nelm:
            is_converged_electronic = True
        if ibrion in [0,1,2,3]:
            if is_job_finished:
                return  is_converged_electronic, is_converged_ionic
        elif ibrion == -1 :
            if is_job_finished:
                return is_converged_electronic, True
        else:
            warnings.warn('Fast convergence check from OUTCAR not yet implemented for IBRION > 3')

    except FileNotFoundError:
        warnings.warn(f'{file} not found')
        return None,None


def get_parameters_from_outcar(tags,file='OUTCAR'):
    paramater_strings = {
        'IBRION': 'ionic relax: 0-MD 1-quasi-New 2-CG',
        'NSW':'number of steps for IOM',
        'NELM':'# of ELM steps'
        }
    
    parameters = {}
    try:
        with open(file,'r') as f:
            lines = f.readlines()
        for line in lines:
            for key in tags:
                target = paramater_strings[key]
                if target in line:
                    characters_without_spaces = [char for char in line.split(' ') if char != '']
                    tag = characters_without_spaces[0]
                    if key != tag:
                        raise ValueError('Target tag and tag read from OUTCAR differ')
                    value = characters_without_spaces[2].replace(';','')
                    parameters[tag] = value
        return parameters
    except FileNotFoundError:
        return None
    
                                

