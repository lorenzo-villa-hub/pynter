#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 16:44:38 2023

@author: villa
"""
import matplotlib

matplotlib.use('Agg') # no graphical output

from pynter.defects.thermodynamics import ThermoData
from pynter.defects.plotter import PressurePlotter

from pynter.testing.core import PynterTest


class TestPressurePlotter(PynterTest):

    # make sure the methods are callable and do not fail
    def test_pressure_plotter(self):
        thermo = ThermoData.from_json(self.get_testfile_path('thermodata_NN_1300K.json'))
        plotter = PressurePlotter()
        plotter.plot_concentrations(thermo)
        plotter.plot_concentrations(thermo,output='all',charges=[0,1])
        plotter.plot_fermi_level(thermo.partial_pressures,thermo.fermi_levels,band_gap=2.48)    