#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 16:44:38 2023

@author: villa
"""

import os
import os.path as op
import unittest
import matplotlib

matplotlib.use('Agg') # no graphical output

from pynter.defects.thermodynamics import ThermoData
from pynter.defects.plotter import PressurePlotter

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/defects/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

# make sure the methods are callable and do not fail
def test_pressure_plotter():
    thermo = ThermoData.from_json(get_path('thermodata_NN_1300K.json'))
    plotter = PressurePlotter()
    plotter.plot_concentrations(thermo)
    plotter.plot_concentrations(thermo,output='all',charges=[0,1])
    plotter.plot_fermi_level(thermo.partial_pressures,thermo.fermi_levels,band_gap=2.48)    