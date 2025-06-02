#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:12:05 2020

@author: lorenzo
"""

import numpy as np
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition


from .core import Chempots, chempot_ideal_gas
from .phase_diagram import PDHandler
from .reservoirs import PressureReservoirs


def get_oxygen_chempot_standard_finite_temperature(temperature,muO_reference=None):
    """
    Get value of oxygen delta mu standard (mu_0(T,Po)) at a speficic temperature
    The data is taken from the following work: 
        Reuter and Scheffler, “Composition, Structure, and Stability of RuO 2 ( 110 ) as a Function of Oxygen Pressure.”
    The value at a specific temperature is extracted from a linear fitting of the behaviour of mu_O with Temperature.
    
    Parameters
    ----------
    temperature : (float)
        Temperature in Kelvin.
    muO_ref : (float)
        Oxygen reference chemical potential (O2 molecule, T = 0K)

    Returns
    -------
    Chemical potential at standard p0 at given temperature.

    """
    T = np.array([100,200,300,400,500,600,700,800,900,1000])
    mu_0 = np.array([-0.08, -0.17,-0.27,-0.38,-0.50,-0.61,-0.73,-0.85,-0.98,-1.10])
    
    coef = np.polyfit(T,mu_0,1)
    poly1d_fn = np.poly1d(coef)
    muO = poly1d_fn(temperature) # delta value
    if muO_reference:
        muO += muO_reference
    return muO


def get_oxygen_chempot_from_pO2(temperature=300,partial_pressure=0.2,muO_reference=None):
    """
    Get oxygen chemical potential (delta) from temperature and partial pressure

    Parameters
    ----------
    temperature : (float), optional
        Temperature in Kelvin.
    partial_pressure : (float)
        Partial pressure.
    muO_ref : (float)
        Oxygen reference chemical potential (O2 molecule, T = 0K)
    Returns
    -------
    (float)
        Value of oxygen chemical potential (delta) at given T and p/p0.
    """
    muO = get_oxygen_chempot_standard_finite_temperature(temperature,muO_reference=muO_reference)
    return chempot_ideal_gas(muO,temperature,partial_pressure)  



def get_pressure_reservoirs_from_phase_diagram(phase_diagram,target_composition,temperature,
                                               extrinsic_chempots_range=None,pressure_range=(-20,10),
                                               interpolation_function=None,
                                               npoints=50,get_pressures_as_strings=False):
    """
    Generate Reservoirs object with a set of different chemical potentials starting from a range of oxygen partial pressure.
    The code distinguishes between 2-component and 3-component phase diagrams.
    
    In the case of a 2-comp PD the other chemical potential is calculated directly from the formation energy of the phase
    with target composition.
    In the case of a 3-comp PD the set of chemical potentials is obtained first calculating the boundary phases starting
    from the constant value of muO, then the arithmetic average between this two points in the stability diagram is taken.
    
    In the case where there are some extrinsic elements not belonging to the PD, they must be added in the 
    extrinsic_chempots_range dictionary ({element:(O_poor_chempot,O_rich_chempot)}). The chempots of extrinsic elements 
    are interpolated linearly between O_poor and O_rich values. 
    
    Parameters
    ----------
    phase_diagram : (PhaseDiagram)
        Pymatgen PhaseDiagram object.
    target_comp : (str or Composition)
        Composition of target phase.
    temperature : (float)
        Temperature in Kelvin.
    extrinsic_chempots_range : (dict)
        Dictionary with chemical potentials of elements not belonging to the PD ({element:(O_poor_chempot,O_rich_chempot)}). 
        The default is None.
    pressure_range : (tuple), optional
        Exponential range in which to evaluate the partial pressure . The default is from 1e-20 to 1e10.
    interpolation_function : (function), optional
        Function to determine the chemical potential for the other elements (the ones that are not oxygen).
        interpolation_function(element,boundary_reservoir):
            The function inputs are the target element and the boundary_reservoirs ({'<label>':{element:value}})
        If None the mean is used. The default is None.    
    npoints : (int), optional
        Number of data points to interpolate the partial pressure with. The default is 50.
    get_pressures_as_strings : (bool), optional
        Get pressure values (keys in the Reservoirs dict) as strings. The default is set to floats.

    Returns
    -------
    reservoirs
        PressureReservoirs object. The dictionary is organized as {partial_pressure:chempots}.
    """
    chempots_dict = {}
    reservoirs = {}
    pd = phase_diagram
    if type(target_composition) == str:
        target_comp = Composition(target_composition)
    else:
        target_comp = target_composition
    partial_pressures = np.logspace(pressure_range[0],pressure_range[1],num=npoints,base=10)
    mu_standard = get_oxygen_chempot_standard_finite_temperature(temperature)
    
    i = 0
    for p in partial_pressures:
        i += 1
        muO = chempot_ideal_gas(mu_standard,temperature=temperature,partial_pressure=p) #delta value
        fixed_chempot = Chempots({'O':muO})
        
        if len(pd.elements) == 2: # 2-component PD case
            pdh = PDHandler(pd)
            mu = pdh.calculate_single_chempot(target_comp, fixed_chempot)
            for el in target_comp.elements:
                if el is not Element('O'):
                    chempots = fixed_chempot.copy()
                    chempots.update({el.symbol:mu})
                
        elif len(pd.elements) == 3: # 3-component PD case
            pdh = PDHandler(pd)
            boundary_res = pdh.get_phase_boundaries_chempots(target_comp, fixed_chempot)
            for el in list(boundary_res.values())[0]:
                if interpolation_function:
                    chempots_dict[el] = interpolation_function(el,boundary_res)
                else:
                    chempots_dict[el] = np.mean(np.array([mu[el] for mu in boundary_res.values()]))
                
        elif len(pd.elements) > 3: 
            raise NotImplementedError('Not implemented for PD with more than 3 components')
        
        chempots = Chempots(chempots_dict)
        chempots_abs = chempots.get_absolute(pdh.mu_refs)
        
        if extrinsic_chempots_range:
            for el in extrinsic_chempots_range:
                mu_el_O_poor,mu_el_O_rich = extrinsic_chempots_range[el][0], extrinsic_chempots_range[el][1]
                chempots_abs[el] = mu_el_O_poor + ((mu_el_O_rich - mu_el_O_poor)/npoints)*i 
                    
        if get_pressures_as_strings:
            p = "%.1g" % p
            p = str(p)
        else:
            p = float("{:.3e}".format(p))
        reservoirs[p] = chempots_abs

        
    return PressureReservoirs(reservoirs,temperature,phase_diagram=pd,are_chempots_delta=False)
                                

def get_oxygen_pressure_reservoirs(oxygen_ref,temperature,pressure_range=(-20,10),npoints=50,get_pressures_as_strings=False):
    """
    Get PressureReservoirs object for oxygen starting from the reference value.

    Parameters
    ----------
    oxygen_ref : (float)
        Absolute chempot of oxygen at 0K.
    temperature : (float), optional
        Temperature.
    pressure_range : (tuple), optional
        Exponential range in which to evaluate the partial pressure . The default is from 1e-20 to 1e10.
    npoints : (int), optional
        Number of data points to interpolate the partial pressure with. The default is 50.
    get_pressures_as_strings : (bool), optional
        Get pressure values (keys in the Reservoirs dict) as strings. The default is set to floats.

    Returns
    -------
    PressureReservoirs object

    """
    reservoirs = {}
    mu_refs = Chempots({'O':oxygen_ref})
    partial_pressures = np.logspace(pressure_range[0],pressure_range[1],num=npoints,base=10)
    mu_standard = get_oxygen_chempot_standard_finite_temperature(temperature)
    for p in partial_pressures:
        mu = {}
        muO = chempot_ideal_gas(mu_standard,temperature=temperature,partial_pressure=p)
        mu.update({'O':oxygen_ref + muO})
        if get_pressures_as_strings:
            p = "%.1g" % p
            p = str(p)
        else:
            p = float("{:.3e}".format(p))
        reservoirs[p] = Chempots(mu)
    
    return PressureReservoirs(reservoirs,temperature,phase_diagram=None,
                              mu_refs=mu_refs,are_chempots_delta=False)

    

        
        
        
        
        