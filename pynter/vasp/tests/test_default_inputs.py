#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:07:14 2023

@author: villa
"""

import os
import os.path as op
from pymatgen.io.vasp.inputs import Poscar, Kpoints

from pynter.vasp.default_inputs import DefaultInputs
from pynter.tests.__init__ import get_structure_Si

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/vasp/tests/test_files')

cfgfile = os.path.join(homedir, ".pynter", "vasp.yml")


structure = get_structure_Si()

def test_get_incar_default():
    di = DefaultInputs(structure=structure, cfgfile=cfgfile)
    incar_dict = di.get_incar_default()
    assert incar_dict["SYSTEM"] == "Si"
    assert incar_dict["ISYM"] == 2

def test_get_kpoints_default():    
    di = DefaultInputs(structure)
    kpoints = di.get_kpoints_default(kppa=1000)
    assert kpoints.kpts == [[8, 8, 8]]
    assert kpoints.kpts_shift == (0, 0, 0)
    
def test_get_kpoints_bs_default():
    di = DefaultInputs(structure)
    kpoints = di.get_kpoints_bs_default(divisions=10,hybrid_mode=False,kppa=1000)
    kpoints = Kpoints.from_string(str(kpoints)) # avoid pymatgen inconsistencies
    kpoints_file = Kpoints.from_file(op.join(test_files_path,'KPOINTS_bs_Si'))
    assert str(kpoints) == str(kpoints_file)
    
    kpoints = di.get_kpoints_bs_default(divisions=10,hybrid_mode=True,kppa=1000)
    kpoints = Kpoints.from_string(str(kpoints)) # avoid pymatgen inconsistencies
    kpoints_file = Kpoints.from_file(op.join(test_files_path,'KPOINTS_bs_hybrid_Si'))
    assert str(kpoints) == str(kpoints_file)


def test_get_potcar():
    di = DefaultInputs(structure=structure, cfgfile=cfgfile)
    potcar = di.get_potcar(potcar_functional='PBE')
    assert len(potcar) == 1
    assert potcar[0].symbol == "Si"


