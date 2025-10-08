
import numpy as np
from scipy.constants import k, h, pi, m_e
from scipy.optimize import bisect

from pymatgen.electronic_structure.dos import Dos, FermiDos, CompleteDos
from pymatgen.electronic_structure.core import Spin
from pymatgen.core.units import kb



def get_carrier_concentrations(dos, fermi_level, temperature, band_gap=None):
    """
    Get carrier concentrations by integrating density of states or using effective masses.
    
    Code from pymatgen.electronic_structure.dos.FermiDos.get_doping
    The function has been modified to return the absolute values of
    hole and electron concentrations.

    Calculate carrier concentrations at a given
    Fermi level and temperature. A simple Left Riemann sum is used for
    integrating the density of states over energy & equilibrium Fermi-Dirac
    distribution.

    Parameters
    ----------
        dos : (dict or Dos)
            Density of states to integrate. Can be provided as density of states D(E)
            or using effective masses.
            Format for effective masses:
                dict with following keys:
                    - "m_eff_h" : holes effective mass in units of m_e (electron mass)
                    - "m_eff_e" : electrons effective mass in units of m_h          
                    - `band_gap` needs to be provided in args
            Formats for explicit DOS:
                dictionary with following keys:
                    - 'energies' : list or np.array with energy values
                    - 'densitites' : list or np.array with total density values
                    - 'structure' : pymatgen Structure of the material, 
                                    needed for DOS volume and charge normalization.
                or a pymatgen Dos object (Dos, CompleteDos or FermiDos).

        fermi_level : (float)
            The Fermi level relative to the VBM in eV.
        temperature : (float)
            The temperature in Kelvin.
        band_gap : (float)
            The band gap in eV. If None is determined from the DOS.

    Returns:
    --------
        h : (float)
            Absolute value of hole concentration in 1/cm^3
        n : (float)
            Absolute value of electron concentration in 1/cm^3
    """   
    if type(dos) == dict:
        if 'energies' in dos and 'densities' in dos:
            E = dos['energies']
            D = dos['densities']
            structure = dos['structure']
            fdos = _get_fermidos_from_data(E=E,D=D,structure=structure,bandgap=band_gap)

        elif 'm_eff_e' in dos and 'm_eff_h' in dos:
            if not band_gap:
                raise ValueError('Band gap must be provided when computing DOS with effective masses')
            
            m_eff_h = dos['m_eff_h'] * m_e
            m_eff_e = dos['m_eff_e'] * m_e
            T = temperature
            
            exponent_h = np.exp(-(fermi_level)/(kb*T))  # fermi_level referenced to VBM
            h = get_dos_from_effective_mass(m_eff_h,T=T) * exponent_h
            
            exponent_e = np.exp(-(band_gap - fermi_level)/(kb*T))
            n = get_dos_from_effective_mass(m_eff_e,T=T) * exponent_e
            return abs(h), abs(n)

    elif type(dos) in  (Dos,FermiDos,CompleteDos):
        fdos = FermiDos(dos=dos,bandgap=band_gap)
    else:
        raise ValueError('DOS must a dictionary (read function docs) or pymatgen Dos object')
    
    _,fdos_vbm = fdos.get_cbm_vbm()
    fermi_level = fdos_vbm + fermi_level

    cb_integral = np.sum(
        fdos.tdos[max(fdos.idx_mid_gap, fdos.idx_vbm + 1) :]
        * f0(fdos.energies[max(fdos.idx_mid_gap, fdos.idx_vbm + 1) :], fermi_level, temperature)
        * fdos.de[max(fdos.idx_mid_gap, fdos.idx_vbm + 1) :],
        axis=0,
    )
    vb_integral = np.sum(
        fdos.tdos[: min(fdos.idx_mid_gap, fdos.idx_cbm - 1) + 1]
        * f0(-fdos.energies[: min(fdos.idx_mid_gap, fdos.idx_cbm - 1) + 1], -fermi_level, temperature)
        * fdos.de[: min(fdos.idx_mid_gap, fdos.idx_cbm - 1) + 1],
        axis=0,
    )
    h = (vb_integral) / (fdos.volume * fdos.A_to_cm ** 3) 
    n = -1*(cb_integral) / (fdos.volume * fdos.A_to_cm ** 3)
    
    return abs(h), abs(n)


def solve_intrinsic_fermi_level(dos,temperature,band_gap,xtol=1e-05):

    def _get_total_q(ef):
        h,n = get_carrier_concentrations(dos=dos,fermi_level=ef,temperature=temperature,band_gap=band_gap)
        return h - n
    
    return bisect(_get_total_q, -1., band_gap + 1.,xtol=xtol)


def f0(E, fermi, T):
    """Fermi-Dirac distribution function.

    Args:
        E (float): Energy in eV.
        fermi (float): The Fermi level in eV.
        T (float): The temperature in kelvin.

    Returns:
        float: The Fermi-Dirac occupation probability at energy E.
    """
    from scipy.constants import value as _constant
    from scipy.special import expit
    exponent = (E - fermi) / (kb * T)
    return expit(-exponent) 


def get_dos_from_effective_mass(m_eff, T):
    """
    Calculate the effective density of states (N_c or N_v)
    for a non-degenerate semiconductor.

    Parameters
    ----------
    m_eff : float
        Effective mass (in units of kg)
    T : float
        Temperature (in Kelvin)

    Returns
    -------
    N : float
        Effective density of states (in cm^-3)
    """
    from scipy.constants import k, h, pi, m_e
    return 2 * ((2 * pi * m_eff * k * T) / (h**2))**(3/2) *1e-06  # units of cm^-3



def _get_fermidos_from_data(E,D,structure,bandgap=None):
    if type(E) == list:
        E = np.array(E)
    if type(D) == list:
        D = np.array(D)

    dos = Dos(energies=E,densities={Spin.up:D})
    return FermiDos(dos=dos,structure=structure,bandgap=bandgap)