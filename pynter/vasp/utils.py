#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 19:20:55 2023

@author: lorenzo
"""
import numpy as np

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