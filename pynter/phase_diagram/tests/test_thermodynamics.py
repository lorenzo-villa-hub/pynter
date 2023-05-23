#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 11:11:39 2023

@author: villa
"""

import os
import os.path as op
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram

from pynter.tools.utils import get_object_from_json
from pynter.phase_diagram.chempots import Chempots, Reservoirs, PDHandler, PressureReservoirs
from pynter.phase_diagram.thermodynamics import OxygenPressure


homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/phase_diagram/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

comp = Composition('NaNbO3')

mu_refs = {'Na': -1.31, 'Nb': -10.1, 'O': -4.95}
chempots_dict = {'Na':-2.26,'Nb':-6.15,'O':-1.92}
mu = Chempots(chempots_dict)
mu_abs = Chempots({'Na':-3.57,'Nb':-16.25,'O':-6.87})


def test_oxygen_pressure():
    oxpr = OxygenPressure(temperature=1300,partial_pressure=0.2)
    oxygen_standard_chempot = oxpr.oxygen_standard_chempot()
    assert oxygen_standard_chempot == -1.4265454545454541
    assert oxpr.chempot_ideal_gas(oxygen_standard_chempot,1300) == -1.516694354091398
    assert oxpr.get_oxygen_chempot_from_pO2() == -1.516694354091398
    ox_res_dict = {
      1e-20: {'O': -10.88},
      3.162277660168379e-13: {'O': -9.91},
      1e-05: {'O': -8.94},
      316.22776601683796: {'O': -7.97},
      10000000000.0: {'O': -7.01}
      }
    ox_res = oxpr.get_oxygen_pressure_reservoirs(mu_abs['O'],npoints=5)
    assert ox_res == PressureReservoirs(ox_res_dict)
    pres_dict = {
     1e-20: {'Na': -2.53, 'Nb': -11.14, 'O': -8.96},
     3.162277660168379e-13: {'Na': -3.07, 'Nb': -13.51, 'O': -7.99},
     1e-05: {'Na': -3.55, 'Nb': -15.94, 'O': -7.02},
     316.22776601683796: {'Na': -4.03, 'Nb': -18.36, 'O': -6.05},
     10000000000.0: {'Na': -4.51, 'Nb': -20.76, 'O': -5.09}
     }
    pd = get_object_from_json(PhaseDiagram,get_path('PD_Na-Nb-O.json'))
    pres = oxpr.get_pressure_reservoirs_from_pd(pd,comp,temperature=1300,npoints=5)
    assert pres == PressureReservoirs(pres_dict)
    assert pres == PressureReservoirs.from_dict(pres.as_dict())
    