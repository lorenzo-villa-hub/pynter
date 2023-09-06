#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 16:44:38 2023

@author: villa
"""
import matplotlib

matplotlib.use('Agg') # no graphical output

from pynter.defects.thermodynamics import ThermoData
from pynter.defects.plotter import ThermodynamicsPlotter

from pynter.testing.core import PynterTest
from pynter.defects.tests.test_thermodynamics import TestDefectThermodynamics


class TestThermodynamicsPlotter(PynterTest):

    # make sure the methods are callable and do not fail
    def test_pressure_plotter(self):
        test_df_thermo = TestDefectThermodynamics()
        test_df_thermo.setUp()
        thermo = test_df_thermo.dt.get_pO2_thermodata(test_df_thermo.pres,temperature=1300)
        plotter = ThermodynamicsPlotter()
        plotter.plot_pO2_vs_concentrations(thermo)
        plotter.plot_pO2_vs_concentrations(thermo,output='all',charges=[0,1])
        plotter.plot_pO2_vs_fermi_level(thermo.partial_pressures,thermo.fermi_levels,band_gap=2.48)    