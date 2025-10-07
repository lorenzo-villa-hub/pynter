#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:56:05 2023

@author: villa
"""
from scipy import constants
from scipy.constants import physical_constants 
import numpy as np
import pandas as pd

from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.core.units import kb

from pynter.defects.analysis import DefectsAnalysis
from pynter.defects.chempots.core import Chempots
from pynter.defects.thermodynamics import DefectThermodynamics
from pynter.tools.utils import get_object_from_json

from pynter.testing.core import PynterTest
#from pynter.defects.tests.test_analysis import TestDefectsAnalysis


e = constants.e

        
class TestDefectThermodynamics(PynterTest):
    
    @classmethod
    def setUpClass(cls):
        self = cls()
        bulk_volume = 909.0035
        d =[{'name': 'Int_O',
            'charge': -2.0,
            'multiplicity': 1,
            'energy_diff': 4.543,
            'bulk_volume': bulk_volume},
            {'name': 'Int_O',
            'charge': -1.0,
            'multiplicity': 1,
            'energy_diff': 1.7746,
            'bulk_volume': bulk_volume},
            {'name': 'Int_O',
            'charge': 0.0,
            'multiplicity': 1,
            'energy_diff': -0.3334,
            'bulk_volume': bulk_volume},
            {'name': 'Sub_P_on_Si',
            'charge': 0.0,
            'multiplicity': 1,
            'energy_diff': 5.6674,
            'bulk_volume': bulk_volume},
            {'name': 'Sub_P_on_Si',
            'charge': 1.0,
            'multiplicity': 1,
            'energy_diff': -1.5959,
            'bulk_volume': bulk_volume},
            {'name': 'Vac_O',
            'charge': 0.0,
            'multiplicity': 1,
            'energy_diff': 9.8442,
            'bulk_volume': bulk_volume},
            {'name': 'Vac_O',
            'charge': 1.0,
            'multiplicity': 1,
            'energy_diff': 8.3845,
            'bulk_volume': bulk_volume},
            {'name': 'Vac_O',
            'charge': 2.0,
            'multiplicity': 1,
            'energy_diff': 4.0104,
            'bulk_volume': bulk_volume},
            {'name': 'Vac_Si',
            'charge': -5.0,
            'multiplicity': 1,
            'energy_diff': 38.4723,
            'bulk_volume': bulk_volume},
            {'name': 'Vac_Si',
            'charge': -4.0,
            'multiplicity': 1,
            'energy_diff': 31.4767,
            'bulk_volume': bulk_volume},
            {'name': 'Vac_Si',
            'charge': -3.0,
            'multiplicity': 1,
            'energy_diff': 28.4873,
            'bulk_volume': bulk_volume}]
        
        df = pd.DataFrame(d)
        vbm = 1.0394
        band_gap = 5.968500000000001
        da = DefectsAnalysis.from_dataframe(df,vbm=vbm,band_gap=band_gap)
        cls.da = da
        mu_O = -6
        mu_P = -5.4133
        mu_Si = -5.4224
        cls.chempots = Chempots({'Si':mu_Si,'P':mu_P,'O':mu_O},ndecimals=2)
        cls.dos = get_object_from_json(CompleteDos, self.get_testfile_path('SiO2-defects/Bulk-2x2x2-supercell/complete_dos.json'))

    
    def test_brouwer_diagram_from_precursors(self):
        da = self.da.filter_entries(elements=['P'],exclude=True)
        da.plot_brouwer_diagram(
                        bulk_dos=self.dos,
                        temperature=1000,
                        precursors={'SiO2':-23.69},
                        oxygen_ref=-4.95,
                        pressure_range=(1e-20,1e10),
                        npoints=50,
                        xtol=1e-15)
        data = da.thermodata
        conc = data.defect_concentrations[5]

        actual = data.partial_pressures[37]
        desired = 449.8
        self.assert_all_close(actual,desired, rtol=1e-03)

        actual = conc.select_concentrations(name='Vac_O',charge=0)[0]['conc']
        desired = 20033806666.685966
        self.assert_all_close(actual,desired, rtol=1e-03)
        
        conc = data.defect_concentrations[21]
        actual = conc.select_concentrations(name='Int_O',charge=-2)[0]['conc']
        desired = 152050.744715061
        self.assert_all_close(actual, desired, rtol=1e-03)

        actual = data.fermi_levels[12]
        desired = 3.0835355572700505
        self.assert_all_close(actual, desired, rtol=1e-03)

    
    def test_doping_diagram(self):
        da = self.da
        da.plot_doping_diagram(
                        variable_defect_specie='P',
                        concentration_range=(1,1e20),
                        chemical_potentials=self.chempots,
                        bulk_dos=self.dos,
                        temperature=1000,
                        npoints=50,
                        xtol=1e-100)
        
        data = da.thermodata
        self.assertEqual(data.variable_defect_specie, 'P')

        actual = data.variable_concentrations[7]
        desired = 719.6856730011522
        self.assert_all_close(actual,desired, rtol=1e-03)

        actual = data.defect_concentrations[7].elemental['P']
        self.assert_all_close(actual, desired, rtol=1e-03)

        conc = data.defect_concentrations[22]
        actual = conc.select_concentrations(name='Int_O',charge=-1)[0]['conc']
        desired = 175.7802582858784
        self.assert_all_close(actual, desired, rtol=1e-03)

        actual = data.fermi_levels[13]
        desired = 2.7925563006401064
        self.assert_all_close(actual, desired, rtol=1e-03)


### more complex functions with quenched species and external defects to be implemented


        
        
        