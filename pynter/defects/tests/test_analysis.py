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
    
    stable_charges = {
         "Int_Si(mult108)": (1.0, 3.6846594634256036),
         "Int_Si(mult54)": (1.0, 3.2615360171135395),
         "Sub_B_on_Si": (0.0, 0.7012506599999799),
         "Sub_P_on_Si": (1.0, -0.25706887286827396),
         "Vac_Si": (0.0, 3.2633273699999723)
         }
    assert da.stable_charges(chempots) == stable_charges 
    
    formation_energies = {
        'Int_Si(mult108)': [(-1.0, 4.649770703425618),
         (0.0, 3.8291055299999766),
         (1.0, 3.6846594634256036)],
        'Int_Si(mult54)': [(-1.0, 4.673740193780213),
         (0.0, 3.8109430399999678),
         (1.0, 3.2615360171135395)],
        'Sub_B_on_Si': [(-2.0, 1.8952842985270144),
         (-1.0, 0.793132301342261),
         (0.0, 0.7012506599999799)],
        'Sub_P_on_Si': [(0.0, 0.2571857899999852),
         (1.0, -0.25706887286827396),
         (2.0, 0.19056955168492085)],
        'Vac_Si': [(-1.0, 3.698017458184358),
         (0.0, 3.2633273699999723),
         (1.0, 3.3491269708159415)]
         }
    
    assert da.formation_energies(chempots == formation_energies)
    
    