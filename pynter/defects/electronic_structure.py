
import numpy as np

from pymatgen.electronic_structure.dos import Dos, FermiDos, CompleteDos
from pymatgen.electronic_structure.core import Spin



def get_carrier_concentrations(dos, fermi_level, temperature, band_gap=None):
    """
    Code from pymatgen.electronic_structure.dos.FermiDos.get_doping
    The function has been modified to return the absolute values of
    hole and electron concentrations.

    Calculate carrier concentrations at a given
    Fermi level and temperature. A simple Left Riemann sum is used for
    integrating the density of states over energy & equilibrium Fermi-Dirac
    distribution.

    -----------
    Parameters:
        dos : (dict or Dos)
            Density of states to integrate. 
            Can either be a dictionary with following keys:
            - 'energies' : list or np.array with energy values
            - 'densitites' : list or np.array with total density values
            - 'structure' : pymatgen Structure of the material, needed for DOS volume normalization.
        fermi_level : (float)
            The Fermi level relative to the VBM in eV.
        temperature : (float)
            The temperature in Kelvin.
        band_gap : (float)
            The band gap in eV. If None is determined from the DOS.

    Returns:
        h : (float)
            Absolute value of hole concentration in 1/cm^3
        n : (float)
            Absolute value of electron concentration in 1/cm^3
    """   
    if type(dos) == dict:
        E = dos['energies']
        D = dos['densities']
        structure = dos['structure']
        fdos = _get_fermidos_from_data(E=E,D=D,structure=structure,bandgap=band_gap)
    elif type(dos) in  (Dos,FermiDos,CompleteDos):
        fdos = FermiDos(dos=dos,bandgap=band_gap)
    else:
        raise ValueError('DOS must be in dict form (read function docs) or pymatgen Dos object')
    
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
    exponent = (E - fermi) / (_constant("Boltzmann constant in eV/K") * T)
    return expit(-exponent) 



def _get_fermidos_from_data(E,D,structure,bandgap=None):
    if type(E) == list:
        E = np.array(E)
    if type(D) == list:
        D = np.array(D)

    dos = Dos(energies=E,densities={Spin.up:D})
    return FermiDos(dos=dos,structure=structure,bandgap=bandgap)