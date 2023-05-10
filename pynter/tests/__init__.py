#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:13:48 2023

@author: villa
"""

from pymatgen.core.structure import Structure

def get_structure_Si():
    structure = Structure.from_dict(
    {'@module': 'pymatgen.core.structure',
     '@class': 'Structure',
     'charge': 0,
     'lattice': {'matrix': [[3.32548851, 0.0, 1.91997169],
       [1.10849617, 3.13530064, 1.91997169],
       [0.0, 0.0, 3.83994338]],
      'pbc': (True, True, True),
      'a': 3.839943374653261,
      'b': 3.839943378813096,
      'c': 3.83994338,
      'alpha': 59.99999998977525,
      'beta': 59.99999995393976,
      'gamma': 60.00000000512866,
      'volume': 40.036809671145996},
     'sites': [{'species': [{'element': 'Si', 'occu': 1}],
       'abc': [0.875, 0.875, 0.875],
       'xyz': [3.879736595, 2.74338806, 6.719900914999999],
       'label': 'Si',
       'properties': {}},
      {'species': [{'element': 'Si', 'occu': 1}],
       'abc': [0.125, 0.125, 0.125],
       'xyz': [0.554248085, 0.39191258, 0.959985845],
       'label': 'Si',
       'properties': {}}]}
        )
    return structure