#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 17:34:15 2023

@author: villa
"""

import os
import os.path as op
from pymatgen.electronic_structure.dos import FermiDos

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots
from pynter.defects.analysis import DefectsAnalysis

from pynter.defects.tests.test_entries import get_defect_entry
from pynter.defects.tests.compare import CompareEntries

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/defects/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

mu_B = -6.6794
mu_P = -5.4133
mu_Si = -5.4224
chempots = Chempots({'Si':mu_Si,'P':mu_P,'B':mu_B})
vbm = 5.8268
band_gap = 0.5729

#dos = get_object_from_json(FermiDos, get_path('Si_DOS.json')) # dos import still not working

entry = get_defect_entry()


def test_defects_analysis():
    da = DefectsAnalysis.from_json(get_path('DA_Si.json'))
    CompareEntries().compare(da.entries[3], entry)
    