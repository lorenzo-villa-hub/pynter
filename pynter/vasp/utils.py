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
    """
    Get charge of vasp calculation from ComputedEntry object. 
    Subtracts the total valence electrons from potcar to 'NELECT' value
    """
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
    """
    Check convergence of VASP calculation from OUTCAR. The following logic is used:
    1. Check for lines in OUTCAR: 
        - Every line with "reached required accuracy" is a succesful electronic step
        - "stopping structural energy minimisation" indicates completed ionic convergence
        - "General timing and accounting information for this job" indicates that the job
          has gone to completion.
    2. Convergence criteria:
        - Electronic convergence is achieved if number_electronic_steps < 'NELM' parameter
        - if 'IBRION' is -1 (no ionic relaxation) also ionic convergence returns True
        - if 'IBRION' is one of (0,1,2,3) ionic convergence is returned based on OUTCAR line
        - If 'IBRION' > 3 ionic convergence returns None
        
    Parameters
    ----------
    file : (str)
        Path to OUTCAR file.

    Returns
    -------
    is_converged_electronic (bool)
    is_converged_ionic (bool)
    """
    try:
        with open(file, "r") as f:
            lines = f.readlines()      
        parameters = {}
        is_converged_electronic, is_converged_ionic, is_job_finished = False,False,False
        n_electronic_steps = 0
        for line in lines:    
            if "reached required accuracy" in line:
                n_electronic_steps += 1
            if "stopping structural energy minimisation" in line:
                is_converged_ionic = True
            if "General timing and accounting informations for this job" in line:
                is_job_finished = True
            param = _get_parameter_from_outcar_line(tags=['IBRION','NELM'],outcar_line=line)
            if param:
                parameters.update(param)
        
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
            return is_converged_electronic, None

    except FileNotFoundError:
        warnings.warn(f'{file} not found')
        return None,None


def get_parameters_from_outcar(tags,file='OUTCAR'):    
    """
    Get parameters from OUTCAR file based on target line search
    """
    parameters = {}
    try:
        with open(file,'r') as f:
            lines = f.readlines()
        for line in lines:
            param = _get_parameter_from_outcar_line(tags=tags,outcar_line=line)
            if param:
               parameters.update(param)
        return parameters
    except FileNotFoundError:
        return None
    
                                
def _get_parameter_from_outcar_line(tags,outcar_line):
    """
    Get parameter value from a line in OUTCAR file based on target line search.
    Lines are defined in the parameter_strings dict.
    """
    paramater_strings = {
        'IBRION': 'ionic relax: 0-MD 1-quasi-New 2-CG',
        'NSW':'number of steps for IOM',
        'NELM':'# of ELM steps'
        }
    params = {}
    for key in tags:
        target = paramater_strings[key]
        if target in outcar_line:
            characters_without_spaces = [char for char in outcar_line.split(' ') if char != '']
            tag = characters_without_spaces[0]
            if key != tag:
                raise ValueError('Target tag and tag read from OUTCAR differ')
            value = characters_without_spaces[2].replace(';','')
            params[tag] = value
    return params

    