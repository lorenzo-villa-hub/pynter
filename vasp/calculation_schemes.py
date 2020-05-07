#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 14:52:51 2020

@author: villa
"""
import os
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import VaspInput, Incar, Poscar, Kpoints, Potcar
from pynter.vasp.default_inputs import DefaultInputs
from pynter.slurm.job_script import ScriptHandler

class CalculationSchemes:
    """
    Class to generate and write input files for different calculation schemes in VASP
    """
    
    def __init__(self,structure=None,incar_settings=None,kpoints=None,potcar=None,job_settings=None,name=None):
        """
        Parameters
        ----------
        structure : (Pymatgen Structure object), optional
            Pymatgen Structure object. If set to None no input parameters can be generated but only lists with step labels.
        incar_settings : (Dict), optional
            Dictionary with incar flags. The default is None. If None the default settings for PBE functional from the DefaultInputs class are used.
        kpoints : (Pymatgen Kpoints object), optional
            Pymatgen Kpoints object. The default is None. If None the default settings from the DefaultInputs class are used.
        potcar : (Pymatgen kpoints object), optional
            Pymatgen kpoints object. The default is None. If None the default settings from the DefaultInputs class are used.
        job_settings : (Dict), optional
            Dictionary with job settings to create job script, parameters are defined in ScrpitHandler class function. The default is None.\n
            If job_settings is None, the 'name' key will be added, the value is the 'name' parameter if provided, if 'name' parameter is \n
            None the value will be: 'no_name'.
        name : (str), optional
            Name for the system to set up scheme for. The default is None.
        """
        
        self.structure = structure if structure else None
        if structure:
            self.incar_settings = incar_settings if incar_settings else DefaultInputs(self.structure).get_incar_default()
            self.kpoints = kpoints if kpoints else DefaultInputs(self.structure).get_kpoints_default()
            self.potcar = potcar if potcar else DefaultInputs(self.structure).get_potcar()
            self.job_settings = job_settings if job_settings else ({'name':name} if name else {'name':'no_name'})
            self.name = name if name else None
            
            if 'name' not in self.job_settings.keys():
                self.job_settings['name'] = self.name
            if 'name' in self.job_settings.keys() and self.name:
                self.job_settings['name'] = self.name
                

       
    def dielectric_properties_electronic(self,scheme_name=None,get_steps_only=False):
        """
        Set calculation for electronic contribution to the dielectric constant (and also dielectric function).
        Uses 'LOPTICS' method in VASP.
        """
        
        scheme = {}
        scheme_name = scheme_name if scheme_name else 'eps-el'
        steps = ['eps-electronic']
        if get_steps_only:
            return steps

        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()  

        val = {}
        for p in self.potcar:
            val[p.element] = p.nelectrons
        incar_settings['NEDOS'] = 2000
        incar_settings['LOPTICS'] = '.TRUE.'              
        nelect = sum([ val[el]*self.structure.composition.as_dict()[el] for el in self.structure.composition.as_dict()]) #number of electrons
        nbands = int((nelect/2 + len(self.structure.sites)/2)*3) # 3*VASP default NBANDS (NELECT/2 +NIONS/2)
        incar_settings['NBANDS'] = nbands

        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name])
        job_settings['add_automation'] = None
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        return Scheme(scheme)
        

    def dielectric_properties_ionic_lcalceps(self,scheme_name=None,get_steps_only=False):
        """
        Set calculation for ionic contribution to the dielectric constant.
        Uses 'LCALCEPS' method in VASP, combined with 'IBRION=6'. Useful for Hybrid calculations where 'LEPSILON' method does not work.
        """
        
        scheme = {}
        scheme_name = scheme_name if scheme_name else 'eps-ion-lcal'
        steps = ['eps-ionic-lcalceps']
        if get_steps_only:
            return steps

        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()  

        val = {}
        for p in self.potcar:
            val[p.element] = p.nelectrons
        incar_settings['NEDOS'] = 2000
        incar_settings['LCALCEPS'] = '.TRUE.'              
        incar_settings['IBRION'] = 6
        incar_settings['NSW'] = 100

        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name])
        job_settings['add_automation'] = None
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        return Scheme(scheme)


    def dielectric_properties_ionic_lepsilon(self,scheme_name=None,get_steps_only=False):
        """
        Set calculation for ionic contribution to the dielectric constant.
        Uses 'LEPSILON' method in VASP, combined with 'IBRION=8'. This method does not work with HSE functionals.
        """
        
        scheme = {}
        scheme_name = scheme_name if scheme_name else 'eps-el-leps'
        steps = ['eps-ionic-lepsilon']
        if get_steps_only:
            return steps

        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()  

        val = {}
        for p in self.potcar:
            val[p.element] = p.nelectrons
        incar_settings['NEDOS'] = 2000
        incar_settings['LEPSILON'] = '.TRUE.'              
        incar_settings['IBRION'] = 8
        incar_settings['NSW'] = 100

        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name])
        job_settings['add_automation'] = None
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        return Scheme(scheme)

    
    def hse_rel(self,scheme_name=None,get_steps_only=False):
        """
        Generates calculation scheme for ionic relaxation for HSE. Steps: \n
            '1-PBE-SCF': Electronic SCF with PBE \n
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE \n
            '3-HSE-SCF': Electronic SCF with HSE \n
            '4-HSE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE
        """
        
        scheme = {}
        scheme_name = scheme_name if scheme_name else 'HSE-rel'
        steps = ['1-PBE-SCF','2-PBE-OPT','3-HSE-SCF','4-HSE-OPT']
        if get_steps_only:
            return steps
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        #set step 1
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'1'])
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'2'])
        
        scheme[steps[1]] = (vaspinput,job_settings)
        
        # set step 3
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'3'])
        
        scheme[steps[2]] = (vaspinput,job_settings)
        
        # set step 4
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'4'])
        
        scheme[steps[3]] = (vaspinput,job_settings)
    
        return Scheme(scheme)
    
    
    def hse_rel_gamma(self,scheme_name=None,get_steps_only=False):
        """
        Generates calculation scheme for structure relaxations for HSE with more intermediate steps. Steps: \n
            '1-PBE-SCF-Gamma': Electronic SCF with PBE only in Gamma point \n
            '2-PBE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE in Gamma point \n
            '3-HSE-SCF-Gamma': Electronic SCF with HSE only in Gamma point \n
            '4-HSE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE in Gamma point \n
            '5-PBE-SCF': Electronic SCF with PBE \n
            '6-HSE-SCF': Electronic SCF with HSE \n
            '7-HSE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE'''
        """
        
        scheme = {}
        scheme_name = scheme_name if scheme_name else'HSE-rel-gamma'
        steps = ['1-PBE-SCF-Gamma','2-PBE-OPT-Gamma','3-HSE-SCF-Gamma',
                  '4-HSE-OPT-Gamma','5-PBE-SCF','6-HSE-SCF','7-HSE-OPT']
        if get_steps_only:
            return steps
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        #set step 1
        incar_settings['NSW'] = 0
        incar_settings['LHFCALC'] = '.FALSE.'
        incar_settings['ISYM'] = 2
        incar = Incar(incar_settings)
        kpoints = Kpoints().gamma_automatic(kpts=(1,1,1))
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'1'])
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'2'])
        
        scheme[steps[1]] = (vaspinput,job_settings)
        
        # set step 3
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'3'])
        
        scheme[steps[2]] = (vaspinput,job_settings)
        
        # set step 4
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'4'])
        
        scheme[steps[3]] = (vaspinput,job_settings)
    
        # set step 5
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'5'])
        
        scheme[steps[4]] = (vaspinput,job_settings)
        
        # set step 6
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'6'])
        
        scheme[steps[5]] = (vaspinput,job_settings)
        
        # set step 7
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'7'])
        
        scheme[steps[6]] = (vaspinput,job_settings) 
        
        return Scheme(scheme)
    
    
    def hse_rel_gamma_extended(self,scheme_name=None,get_steps_only=False):
        """
        Generates calculation scheme for structure relaxations for HSE with more intermediate steps. Steps: \n
            '1-PBE-SCF-Gamma': Electronic SCF with PBE only in Gamma point \n
            '2-PBE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE in Gamma point \n
            '3-HSE-SCF-Gamma': Electronic SCF with HSE only in Gamma point \n
            '4-HSE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE in Gamma point \n
            '5-PBE-SCF': Electronic SCF with PBE \n
            '6-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE \n
            '7-HSE-SCF': Electronic SCF with HSE \n
            '8-HSE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE'''
        """
        scheme = {}
        scheme_name = scheme_name if scheme_name else'HSE-rel-gamma-ext'
        steps = ['1-PBE-SCF-Gamma','2-PBE-OPT-Gamma','3-HSE-SCF-Gamma',
                  '4-HSE-OPT-Gamma','5-PBE-SCF','6-PBE-OPT','7-HSE-SCF','8-HSE-OPT']
        if get_steps_only:
            return steps
        
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        #set step 1
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = Kpoints().gamma_automatic(kpts=(1,1,1))
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'1'])
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'2'])
        
        scheme[steps[1]] = (vaspinput,job_settings)
        
        # set step 3
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'3'])
        
        scheme[steps[2]] = (vaspinput,job_settings)
        
        # set step 4
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'4'])
        
        scheme[steps[3]] = (vaspinput,job_settings)
    
        # set step 5
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'5'])
        
        scheme[steps[4]] = (vaspinput,job_settings)

        # set step 6
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'6'])
        
        scheme[steps[5]] = (vaspinput,job_settings)
        
        # set step 7
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'7'])
        
        scheme[steps[6]] = (vaspinput,job_settings)
        
        # set step 8
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'8'])
        
        scheme[steps[7]] = (vaspinput,job_settings) 
        
        return Scheme(scheme)
                      

    def hse_rel_short(self,scheme_name=None,get_steps_only=False):
        """
        Generates calculation scheme for ionic relaxation for HSE (short version i.e. no PBE preliminary calculation). Steps: \n
            '1-PBE-SCF': Electronic SCF with HSE \n
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE
        """

        scheme = {}
        scheme_name = scheme_name if scheme_name else 'HSE-rel-short'
        steps = ['1-HSE-SCF','2-HSE-OPT']        
        if get_steps_only:
            return steps
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        # set step 1
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar        
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['nodes'] = 8
        job_settings['timelimit'] = '72:00:00'
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'1'])
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['nodes'] = 8
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'2'])
        
        scheme[steps[1]] = (vaspinput,job_settings)
    
        return Scheme(scheme)  

    
    def hse_vol_rel(self,scheme_name=None,get_steps_only=False):
        """
        Generates calculation scheme for structure relaxations for HSE including cell volume relaxation. Steps: \n
            '1-PBE-SCF': Electronic SCF with PBE
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE
            '3-HSE-SCF': Electronic SCF with HSE
            '4-HSE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE
            '5-HSE-VOPT': Cell volume relaxation and ionic relaxation with HSE (ISIF=3)
        """
        
        scheme = {}
        scheme_name = scheme_name if scheme_name else 'HSE-Vrel'
        steps = ['1-PBE-SCF','2-PBE-OPT','3-HSE-SCF','4-HSE-OPT','5-HSE-VOPT']
        if get_steps_only:
            return steps
        
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        #set step 1
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'1'])
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'2'])
        
        scheme[steps[1]] = (vaspinput,job_settings)
        
        # set step 3
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'3'])
        
        scheme[steps[2]] = (vaspinput,job_settings)
        
        # set step 4
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'4'])
        
        scheme[steps[3]] = (vaspinput,job_settings)
        
        # set volume relaxation
        job_settings = self.job_settings.copy()
        incar_settings['ISIF'] = 3
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'5'])
        
        scheme[steps[4]] = (vaspinput,job_settings)

        return Scheme(scheme)
    
    
    def hubbard_tuning(self,specie,ldauu_dict=None,u_range=(1,10),scheme_name=None):
        """
        Generates Scheme for many calculations using different Hubbard corrections

        Parameters
        ----------
        specie : (str)
            Symbol of atomic specie to apply changes on U parameter.
        ldauu_dict : (dict), optional
            Dictionary of U parameters. Keys are pymatgen Element objects, values are the associated U parameters. The default is None.
        u_range : (Tuple), optional
            Tuple with range of U parameter to scan. The default is (1,10).
        scheme_name : (str), optional
            Name for Scheme, if None the name is set to "U_tuning". The default is None.

        Returns
        -------
        Scheme object
        """
        
        scheme = {}
        scheme_name = scheme_name if scheme_name else 'U_tuning'
        
        if ldauu_dict is None:
            ldauu_dict={}
            for el in self.structure.composition.elements:
                ldauu_dict[el] = 0
        else:
            if len(ldauu_dict) != len(self.structure.composition.elements):
                raise ValueError('size of "ldauu_dict" needs to be the same as the number of species in Structure')
                
        for u in range(u_range[0],u_range[1]+1):
            
            incar_settings = self.incar_settings.copy()
            job_settings = self.job_settings.copy()
            
            ldauu_dict[Element(specie)] = u
            incar_settings['LDAUU'] = ' '.join([str(ldauu_dict[el]) for el in ldauu_dict])
            
            incar = Incar(incar_settings)
            kpoints = self.kpoints
            poscar = Poscar(self.structure)
            potcar = self.potcar
            vaspinput = VaspInput(incar,kpoints,poscar,potcar)
            job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,f'U_{u}'])
            
            scheme[f'U_{u}'] = (vaspinput,job_settings)
            
        return Scheme(scheme)
    
    
    def pbe_electronic_structure(self,kmesh_dos=3, kpoints_bs=None,scheme_name=None,get_steps_only=False):
        """
        Generates calculation scheme for electronic structure calculations for PBE (in general non hybrid functionals) Steps: \n
            '1-PBE-SCF': Electronic SC \n
            '2-PBE-DOS': DOS calculation with 2000 number of gridpoints on which the DOS is evaluated (NEDOS=2000) and increased k-mesh, ISMEAR is set to -5.
            '3-PBE-BS': Bandstructure calculation
        Parameters
        ----------
        kmesh_dos : (Int), optional
            Multiplier for coefficients of the k-mesh for the DOS with respect to first step SCF calculation. The default is 3.
        kpoints_bs : (Pymatgen Kpoints object) , optional
            Pymatgen Kpoints object for the path in BS calculation. The default is None. If None the default high symmetry path \n
            from Pymatgen class HighSymmKpath is obtained from the input Structure with 10 points between high symm k-points. '
            
        Returns
        -------
        Scheme object.
        """
     
        scheme = {}
        scheme_name = scheme_name if scheme_name else 'PBE-el-str'
        steps = ['1-PBE-SCF','2-PBE-DOS','3-PBE-BS']
        if get_steps_only:
            return steps
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        incar_settings['NSW'] = 0
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'1'])
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        # set DOS calculation
        job_settings = self.job_settings.copy()
        incar_settings['ISTART'] = 1
        incar_settings['ICHARG'] = 1
        incar_settings['NEDOS'] = 2000
        incar_settings['ISMEAR'] = -5
        incar = Incar(incar_settings)
        
        # multiply by 3 coeff of k-mesh
        style, kpts, kpts_shift = kpoints.style, kpoints.kpts, kpoints.kpts_shift
        kpts_new = []
        kpts_new.append([k*kmesh_dos for k in kpts[0]])
        kpoints = Kpoints(style=style,kpts=kpts_new,kpts_shift=kpts_shift)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'2'])
        
        scheme[steps[1]] = (vaspinput,job_settings)

        # set BS calculation
        job_settings = self.job_settings.copy()
        incar_settings['ICHARG'] = 11
        incar_settings['LORBIT'] = 11
        incar_settings['ISMEAR'] = 0
        if kpoints_bs:
            kpoints = kpoints_bs
        else:
            kpoints = Kpoints().automatic_linemode(10,HighSymmKpath(self.structure))
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'3'])
        
        scheme[steps[2]] = (vaspinput,job_settings)  
        
        return Scheme(scheme)
                 
    
    def pbe_rel(self,scheme_name=None,get_steps_only=False):
        """
        Generates calculation scheme for ionic relaxation with PBE functional. Steps: \n
            '1-PBE-SCF': Electronic SCF \n
            '2-PBE-OPT': Ionic relaxation with maximum 100 ionic steps (NSW=100), Energy change between electronic steps of 1e-05 (EDIFF=1e-05), \n
                         Force convergence criterion of 0.05 eV/Amstrong (EDIFFG=-0.05).

        Returns
        -------
        Scheme object
        """        

        scheme = {}
        scheme_name = scheme_name if scheme_name else 'PBE-rel'
        steps = ['1-PBE-SCF','2-PBE-OPT']
        if get_steps_only:
            return steps
        
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        incar_settings['NSW'] = 0
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'1'])
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        # set ionic relaxation
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'2'])
        
        scheme[steps[1]] = (vaspinput,job_settings)
        
        return Scheme(scheme)
    
    
    def pbe_vol_rel(self,scheme_name=None,get_steps_only=False):
        """
        Generates calculation scheme for ionic and cell relaxation with PBE functional. Steps: \n
            '1-PBE-SCF': Electronic SC \n
            '2-PBE-OPT': Ionic relaxation with maximum 100 ionic steps (NSW=100), Energy change between electronic steps of 1e-05 (EDIFF=1e-05), \n
                         Force convergence criterion of 0.05 eV/Amstrong (EDIFFG=-0.05) \n
            '3-PBE-VOPT': Same parameters for ionic relaxation and relaxation of the cell (ISIF=3).

        Returns
        -------
        Scheme object
        """
        
        scheme = {}
        scheme_name = scheme_name if scheme_name else 'PBE-Vrel'
        steps = ['1-PBE-SCF','2-PBE-OPT','3-PBE-VOPT']
        if get_steps_only:
            return steps
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        incar_settings['NSW'] = 0
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'1'])
        
        scheme[steps[0]] = (vaspinput,job_settings)
        
        # set ionic relaxation
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'2'])
        
        scheme[steps[1]] = (vaspinput,job_settings)
        
        # set volume relaxation
        job_settings = self.job_settings.copy()
        incar_settings['ISIF'] = 3
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,'3'])
        
        scheme[steps[2]] = (vaspinput,job_settings)
        
        return Scheme(scheme)
                 
        
        
class Scheme:
    """ 
    Class to organize data for a calculation scheme. It is unlickely to generate this class manually.
    """
    
    def __init__(self,scheme):
        """
        Parameters
        ----------
        scheme : (Dict)
            Dictionary with step names as keys, and tuple with VaspInput objects (pymatgen.io.vasp.inputs) and dictionary with job settings as values. \n
            ( { step_name: (VaspInput object, job_settings dict) } ).
        """
        self.scheme = scheme
    
    @property
    def steps(self):
        """
        List with steps in the calculation scheme.

        Returns
        -------
        steps : List
            List of steps.
        """
        steps = [s for s in self.scheme]
        return steps
    
    def get_job_settings(self,step):
        """
        Get dictionary with job settings of a given step.

        Parameters
        ----------
        step : (str)
            Step to get associated job settings.

        Returns
        -------
        Dict
            dictionary of job settings.

        """
        script_handler = ScriptHandler(**self.scheme[step][1])
        return script_handler.settings

        
    def get_parameters(self,step):
        """
        Get tuple with Pymatgen VaspInput object and job settings dictionary of a given step.

        Parameters
        ----------
        step : (str)
            Step to get associated parameters.

        Returns
        -------
        Tuple
            (VaspInput object, job settings Dictionary).

        """
        return self.scheme[step]

    
    def get_vaspinput(self,step):
        """
        Get Pymatgen VaspInput object of a given step.

        Parameters
        ----------
        step : (str)
            Step to get associated VaspInput.

        Returns
        -------
        Pymatgen VaspInput object (pymatgen.io.vasp.inputs)

        """
        return self.scheme[step][0]


    def set_job_settings(self,step,job_settings):
        """
        Set job settings of a given step

        Parameters
        ----------
        step : (str)
            Step name.
        job_settings : (Dict)
            Dictionary with job settings.
        """
        self.scheme[step] = (self.get_vaspinput(step),job_settings)
        
    def set_parameters(self,step,vaspinput,job_settings):
        """
        Set both VaspInput and job settings of a given step

        Parameters
        ----------
        step : (str)
            Step name.
        vaspinput: (Pymatgen VaspInput object)
        
        job_settings : (Dict)
            Dictionary with job settings.
        """        
        self.scheme[step] = (vaspinput,job_settings)

    def set_vaspinput(self,step,vaspinput):
        """
        Set VaspInput of a given step

        Parameters
        ----------
        step : (str)
            Step name.
        vaspinput: (Pymatgen VaspInput object)
        """         
        self.scheme[step] = (vaspinput,self.get_job_settings(step))
        

    def write_scheme(self, path=None, make_dir_if_not_present=True):
        """
        Function to write input files from previously generated Scheme object.

        Parameters
        ----------
        scheme : (Scheme object)
            Previously generated Scheme object. It is unlikely that you generate this object manually.
        path : (str), optional
            Path to save calculation inputs. The default is None. If None files are saved in work dir.
        make_dir_if_not_present : (Bool), optional
            Create directory from given path is not present. The default is True.
        """
        
        for step in self.steps:
            complete_path = os.path.join(path, step) if path else step
            vaspinput=self.get_vaspinput(step)
            vaspinput.write_input(complete_path,make_dir_if_not_present=make_dir_if_not_present)
            script_handler = ScriptHandler(**self.get_job_settings(step))
            script_handler.write_script(path=complete_path)
            
        return
