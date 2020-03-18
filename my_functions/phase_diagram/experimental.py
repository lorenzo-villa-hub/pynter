#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:12:05 2020

@author: lorenzo
"""

import numpy as np

class ChempotExperimental:
    
    def __init__(self,temperature=300,partial_pressure=1):
        """
        Class to extract and organize data regarding experimental chemical potentials

        Parameters
        ----------
        temperature : (float), optional
            Temperature. The default is 300.
        partial_pressure : (float), optional
            Partial pressure (p/p0). The default is 1.
        """
        self.temperature = temperature
        self.partial_pressure = partial_pressure


    def chempot_from_ideal_gas(self, mu0, temperature=None,partial_pressure=None):
        """
        Get chemical potential at a given temperature and partial pressure. The chemical potential in standard conditions (mu0)
        has to be know.

        Parameters
        ----------
        mu0 : (float)
            Value of chemical potential in standard conditions.
        temperature : (float), optional
            Temperature. The default is None. If None the value initialized with the class is used.
        partial_pressure : (float), optional
            Partial pressure. The default is None. If None the value initialized with the class is used.

        Returns
        -------
        chempot : (float)
            Value of chemical potential at given T and p/p0.
        """
        temperature = temperature if temperature else self.temperature
        partial_pressure = partial_pressure if partial_pressure else self.partial_pressure
        kb = 8.6173324e-5  # eV / K
        chempot = mu0 + 0.5*kb*temperature*np.log(partial_pressure)
        return chempot
        

        
    def oxygen_standard_chempot(self,temperature=None):
        """
        Get value of oxygen delta mu standard (mu_0(T,Po)) at a speficic temperature
        The data is taken from the following work: 
            Reuter and Scheffler, “Composition, Structure, and Stability of RuO 2 ( 110 ) as a Function of Oxygen Pressure.”
        The value at a specific temperature is extracted from a linear fitting of the behaviour of mu_O with Temperature.
        
        Parameters
        ----------
        temperature : (float)

        Returns
        -------
        Chemical potential at standard p0 at given temperature.

        """
        temperature = temperature if temperature else self.temperature
        T = np.array([100,200,300,400,500,600,700,800,900,1000])
        mu_0 = np.array([-0.08, -0.17,-0.27,-0.38,-0.50,-0.61,-0.73,-0.85,-0.98,-1.10])
        
        coef = np.polyfit(T,mu_0,1)
        poly1d_fn = np.poly1d(coef)
        return poly1d_fn(temperature)