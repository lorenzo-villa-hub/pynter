#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:59:10 2020

@author: villa
"""

from pymatgen.io.vasp.inputs import Kpoints, Poscar, Potcar, VaspInput, Incar
from pynter.vasp.__init__ import load_vasp_default
import os

homedir = os.getenv("HOME")
cfgfile = os.path.join(homedir,'.pynter','vasp.yml')


class DefaultInputs:
    """
    Class to generate VASP input files with default parameters. Use with extreme care.    
    """
    
    def __init__(self,structure=None,cfgfile=cfgfile):
        """
        Parameters
        ----------
        structure : Pymatgen Structure object. The default is None. If no Structure is given only default INCAR can be generated.
        cfgfile : File with default VASP parameters. The default is "~./pynter.vasp.yml".
        """
        self._structure = structure if structure else None 
        defaults = load_vasp_default(cfgfile=cfgfile)
        
        self.default_potcar_symbols =  defaults['default_potcar_symbols']
        if self.structure:
            self.potcar_symbols = [self.default_potcar_symbols[el.symbol] for el in self.structure.composition.elements]
        self.incar_default_flags = defaults['incar_default']
        
    @property
    def structure(self):
        return self._structure
    
    @structure.setter
    def structure(self,structure):
        self._structure = structure

    def get_incar_default(self, xc='PBE', ldauu=None, aexx=0.25):
        """
        Get dictionary with default INCAR parameters

        Parameters
        ----------
        xc : (str), optional
            Which functional to use('PBE','PBE+U','LDA','LDA+U',HSE06' available). The default is 'PBE'.
        ldauu : (List), optional
            List of integers for U parameter to be given in order as required by VASP. The default is None. If None no 'LDAUU' flag is written.
        aexx : (Int), optional
            Value of Hartree-Fock contribution for HSE calculation. The default is 0.25.

        Returns
        -------
        incar_default_flags : Dict
            Dictionary with default INCAR parameters.

        """
        system = self.structure.composition.reduced_formula if self.structure else 'No system info'
          
        incar_default_flags = self.incar_default_flags              
        incar_default_flags["SYSTEM"] = system                   
                     
        if xc == 'PBE' or xc == 'LDA':
            incar_default_flags.update({
                    f"#### Default {xc}: system" : system , 
                    "ISYM":2})
            
        if xc == 'PBE+U' or xc == 'LDA+U':
            incar_default_flags.update({
                    f"#### Default {xc}: system" : system,                 
                    "ISYM":2,
                    "LDAU" : ".TRUE.",
                    "LDAUTYPE": 2,
                    "LDAUPRINT": 2,
                })
            if ldauu:
                incar_default_flags.update({"LDAUU": ' '.join([str(u) for u in ldauu])}) # if LDAUU is not present VASP will put all zeros
            
        if xc == 'HSE06':
            incar_default_flags.update({
                    "#### Default HSE06: system" : system, 
                    "LHFCALC" : ".TRUE.",
                    "HFSCREEN": 0.2,
                    "PRECFOCK": "Fast",
                    "AEXX": aexx,
                    "ISYM":3
                    })
    
        return incar_default_flags


    def get_vasp_input(self, xc='PBE', ldauu=None, aexx=0.25, kppa=1000,
                       potcar_symbols=None, potcar_functional='PBE'):
        """
        Get pymatgen VaspInput object with complete set of default inputs.

        Parameters
        ----------
        xc : (str), optional
            Which functional to use('PBE','PBE+U','HSE06' available). The default is 'PBE'.
        ldauu : (List), optional
            List of integers for U parameter to be given in order as required by VASP. The default is None. If None no 'LDAUU' flag is written.
        aexx : (Int), optional
            Value of Hartree-Fock contribution for HSE calculation. The default is 0.25.
        kppa : (Int), optional
            Value of k-points density per atom to generate automatic k-mesh as done by Pymatgen. The default is 1000.
        potcar_symbols : List, optional
            List of strings with symbols of POTCARs. The default is None. If None the default symbols from the Materials Project \n
            database are used for every specie in the given Structure object.
        potcar_functional : (str), optional
            Functional to be used, as defined by Pymatgen Potcar class. The default is 'PBE'.

        Returns
        -------
        Pymatgen VaspInput object.

        """
        
        incar = Incar(self.get_incar_default(xc=xc,ldauu=ldauu,aexx=aexx))
        kpoints = self.get_kpoints_default(kppa=kppa)
        poscar = self.get_poscar()
        potcar = self.get_potcar(potcar_symbols=potcar_symbols,potcar_functional=potcar_functional)
        
        return VaspInput(incar,kpoints,poscar,potcar)
        
        
    def get_kpoints_default(self,kppa=1000):
        """
        Get Pymatgen Kpoints object with automatic Gamma centered k-mesh generation from given Structure object.

        Parameters
        ----------
        kppa : (Int), optional
            Value of k-points density per atom to generate automatic k-mesh as done by Pymatgen. The default is 1000.

        Returns
        -------
        Pymatgen Kpoint object.
        """
        
        if self.structure:
            structure = self.structure
        else:
            raise ValueError('Structure object needs to be given to obtain pymatgen automatic kpoint mesh')
         
        return Kpoints().automatic_gamma_density(structure,kppa)
    
    
    def get_poscar(self):
        """
        Get Pymatgen Poscar object from given Structure object.        

        Returns
        -------
        Pymatgen Poscar object
        """
        
        if self.structure:
            pass
        else:
            raise ValueError('Structure object needs to be given to get Poscar onject')
        
        return Poscar(self.structure)
        
     
    def get_potcar(self,potcar_symbols=None,potcar_functional='PBE'):
        """
        Get Pymatgen Potcar object. The PMG_VASP_PSP directory has to be defined in .pmgrc.yaml file as required by Pymatgen.
        Parameters
        ----------
        potcar_symbols : List, optional
            List of strings with symbols of POTCARs. The default is None. If None the default symbols from the Materials Project \n
            database are used for every specie in the given Structure object.
        potcar_functional : (str), optional
            Functional to be used, as defined by Pymatgen Potcar class. The default is 'PBE'.

        Returns
        -------
        Pymatgen Potcar object.
        """
                
        if potcar_symbols:
            potcar = Potcar(symbols=potcar_symbols,functional=potcar_functional)
            self.potcar_symbols = potcar_symbols
        else:
            if self.structure:
                potcar_symbols = self.potcar_symbols
            else:
                raise ValueError('potcar symbols or Structure object need to be given to get Potcar onject')
            potcar = Potcar(symbols=potcar_symbols,functional=potcar_functional)
        
        return potcar
            
