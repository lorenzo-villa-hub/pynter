#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:59:10 2020

@author: villa
"""

from pymatgen.io.vasp.inputs import Kpoints, Poscar, Potcar, VaspInput, Incar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

from pynter import SETTINGS


class DefaultInputs:
    """
    Class to generate VASP input files with default parameters. Use with extreme care.    
    """
    
    def __init__(self,structure=None,sort_structure=True):
        """
        Parameters
        ----------
        structure : Pymatgen Structure object. The default is None. If no Structure is given only default INCAR can be generated.
        sort_structure (bool) : Sort structure to avoid POSCAR/POTCAR inconsistencies. The default is True.
        """
        if structure and sort_structure:
            structure.sort()
        self._structure = structure if structure else None 
        defaults = SETTINGS['vasp']
        
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
                    "ISYM":2})
            
        if xc == 'PBE+U' or xc == 'LDA+U':
            incar_default_flags.update({
                    "ISYM":2,
                    "LDAU" : ".TRUE.",
                    "LDAUTYPE": 2,
                    "LDAUPRINT": 2,
                })
            if ldauu:
                incar_default_flags.update({"LDAUU": ' '.join([str(u) for u in ldauu])}) # if LDAUU is not present VASP will put all zeros
            
        if xc == 'HSE06':
            incar_default_flags.update({
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
    
  
    def get_kpoints_bs_default(self,divisions=10,hybrid_mode=False,kppa=1000):
        """
        Get Pymatgen Kpoints object for band structure calculations.
        The default high symmetry path for PBE from Pymatgen class HighSymmKpath is obtained from 
        the input Structure with 10 points between high symm k-points.
        For Hybrid BS the calculations need to be selfconsistent, the KPOINTS file is a list of single kpoints.
        The irriducible kpoints with weight != 0 are obtained with Pymatgen's SpaceGroupAnalyzer.
        The additional points with weight = 0 on the high symm path are obtained with the HighSymmPath class
        like for the PBE case.'

        Parameters
        ----------
        divisions : int, optional
            Number of points in between high symmetry points. The default is 10.
        hybrid_mode : (bool), optional
            Kpoints for hybrid calculation. The default is False.
        kppa : (Int), optional
            Value of k-points density per atom to generate automatic k-mesh as done by Pymatgen. The default is 1000.

        Returns
        -------
        Pymatgen Kpoint object.
        """
        if hybrid_mode: 
            # code copy/pasted (slightly adjusted) from pymatgen.io.vasp.sets.MPHSEBSSet
            # The whole input set was not needed but only the kpoints generation in the kpoints @property
            kpts = []
            weights = []  
            all_labels = []  
            structure = self.structure
    
            # for both modes, include the Uniform mesh w/standard weights
            grid = Kpoints().automatic_gamma_density(structure, kppa).kpts
            ir_kpts = SpacegroupAnalyzer(structure, symprec=0.1).get_ir_reciprocal_mesh(grid[0])
            for k in ir_kpts:
                kpts.append(k[0])
                weights.append(int(k[1]))
                all_labels.append('')
       
            #add the symmetry lines w/zero weight
            kpath = HighSymmKpath(structure)
            frac_k_points, labels = kpath.get_kpoints(line_density=divisions,coords_are_cartesian=False)

            for k, f in enumerate(frac_k_points):
                kpts.append(f)
                weights.append(0.0)
                all_labels.append(labels[k])
    
            comment = "HSE run along symmetry lines"   
            kpoints =  Kpoints(comment=comment,style=Kpoints.supported_modes.Reciprocal,num_kpts=len(kpts),
                kpts=kpts,kpts_weights=weights,labels=all_labels)

        else:
            kpoints = Kpoints().automatic_linemode(10,HighSymmKpath(self.structure))
    
        return kpoints
    
    
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
            
