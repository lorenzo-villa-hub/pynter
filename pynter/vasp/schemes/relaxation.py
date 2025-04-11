#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 14:52:51 2020

@author: villa
"""

from pymatgen.core.periodic_table import Element

from pynter.vasp.schemes.core import InputSets

  
class RelaxationSchemes(InputSets):
    """
    Class to generate and write input files for different calculation schemes in VASP
    """

    def hse_relaxation(self,scheme_name='HSE-rel'):
        """
        Generates calculation scheme for ionic relaxation for HSE. Steps: \n
            '1-PBE-SCF': Electronic SCF with PBE \n
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE \n
            '3-HSE-SCF': Electronic SCF with HSE \n
            '4-HSE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE
            
        Returns:
            List of VaspJob objects
        """        
        stepnames = ['1-PBE-SCF','2-PBE-OPT','3-HSE-SCF','4-HSE-OPT']
        jobs = []
        
        sn = 1 #set step number 1
        vaspjob = self.pbe_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 2
        vaspjob = self.pbe_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 3 
        vaspjob = self.hse_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4
        vaspjob = self.hse_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)

        return jobs
    
    
    def hse_relaxation_gamma(self,scheme_name='HSE-rel-gamma'):
        """
        Generates calculation scheme for structure relaxations for HSE with more intermediate steps. Steps: \n
            '1-PBE-SCF-Gamma': Electronic SCF with PBE only in Gamma point \n
            '2-PBE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE in Gamma point \n
            '3-HSE-SCF-Gamma': Electronic SCF with HSE only in Gamma point \n
            '4-HSE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE in Gamma point \n
            '5-PBE-SCF': Electronic SCF with PBE \n
            '6-HSE-SCF': Electronic SCF with HSE \n
            '7-HSE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE'''

        Returns:
            List of VaspJob objects
        """
        stepnames = ['1-PBE-SCF-Gamma','2-PBE-OPT-Gamma','3-HSE-SCF-Gamma',
                  '4-HSE-OPT-Gamma','5-PBE-SCF','6-HSE-SCF','7-HSE-OPT']
        jobs = []
        
        sn = 1 #set step number 1
        vaspjob = self.pbe_scf_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 2 
        vaspjob = self.pbe_ionic_rel_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 3 
        vaspjob = self.hse_scf_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4 
        vaspjob = self.hse_ionic_rel_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)

        sn = 5 
        vaspjob = self.pbe_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 6 
        vaspjob = self.hse_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 7
        vaspjob = self.hse_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)

        return jobs
    
    
    def hse_relaxation_gamma_extended(self,scheme_name='HSE-rel-gamma-ext'):
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

        Returns:
            List of VaspJob objects
        """
        stepnames = ['1-PBE-SCF-Gamma','2-PBE-OPT-Gamma','3-HSE-SCF-Gamma',
                  '4-HSE-OPT-Gamma','5-PBE-SCF','6-PBE-OPT','7-HSE-SCF','8-HSE-OPT']
        jobs = []
        
        sn = 1 #set step number 1
        vaspjob = self.pbe_scf_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 2
        vaspjob = self.pbe_ionic_rel_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 3
        vaspjob = self.hse_scf_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4
        vaspjob = self.hse_ionic_rel_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)

        sn = 5 
        vaspjob = self.pbe_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 6 
        vaspjob = self.pbe_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 7 
        vaspjob = self.hse_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 8 
        vaspjob = self.hse_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)

        return jobs
                      

    def hse_relaxation_short(self,scheme_name='HSE-rel-short'):
        """
        Generates calculation scheme for ionic relaxation for HSE (short version i.e. no PBE preliminary calculation). Steps: \n
            '1-PBE-SCF': Electronic SCF with HSE \n
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE

        Returns:
            List of VaspJob objects
        """
        stepnames = ['1-HSE-SCF','2-HSE-OPT']        
        jobs = []
        
        sn = 1 
        vaspjob = self.hse_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 2
        vaspjob = self.hse_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        return jobs

    
    def hse_vol_relaxation(self,scheme_name='HSE-Vrel'):
        """
        Generates calculation scheme for structure relaxations for HSE including cell volume relaxation. Steps: \n
            '1-PBE-SCF': Electronic SCF with PBE
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE
            '3-HSE-SCF': Electronic SCF with HSE
            '4-HSE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE
            '5-HSE-VOPT': Cell volume relaxation and ionic relaxation with HSE (ISIF=3)

        Returns:
            List of VaspJob objects
        """
        stepnames = ['1-PBE-SCF','2-PBE-OPT','3-HSE-SCF','4-HSE-OPT','5-HSE-VOPT']
        jobs = []
        
        sn = 1 
        vaspjob = self.pbe_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 2
        vaspjob = self.pbe_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 3 
        vaspjob = self.hse_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4 
        vaspjob = self.hse_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)

        sn = 5
        vaspjob = self.hse_ionic_vol_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        return jobs
    
    
    def hubbard_tuning(self,specie,ldauu_dict=None,u_range=(1,10),scheme_name='U_tuning'):
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
            Name for scheme. The default is "U_tuning".

        Returns
        -------
        List of VaspJob objects
        """
        jobs = []
        
        if ldauu_dict is None:
            ldauu_dict={}
            for el in self.structure.composition.elements:
                ldauu_dict[el] = 0
        else:
            if len(ldauu_dict) != len(self.structure.composition.elements):
                raise ValueError('size of "ldauu_dict" needs to be the same as the number of species in Structure')
                
        for u in range(u_range[0],u_range[1]+1):
                        
            ldauu_dict[Element(specie)] = u
            stepname = f'U{u}'
            vaspjob = self.hubbard(ldauu_dict,add_to_job_name=stepname,add_to_job_path=stepname)
            jobs.append(vaspjob)
            
        return jobs

    
    def pbe_relaxation(self,scheme_name='PBE-rel'):
        """
        Generates calculation scheme for ionic relaxation with PBE functional. Steps: \n
            '1-PBE-SCF': Electronic SCF \n
            '2-PBE-OPT': Ionic relaxation with maximum 100 ionic steps (NSW=100), Energy change between electronic steps of 1e-05 (EDIFF=1e-05), \n
                         Force convergence criterion of 0.05 eV/Amstrong (EDIFFG=-0.05).

        Returns
        -------
        List of VaspJob objects
        """        
        stepnames = ['1-PBE-SCF','2-PBE-OPT']
        jobs = []
              
        sn = 1
        vaspjob = self.pbe_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 2 
        vaspjob = self.pbe_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        return jobs


    def pbe_relaxation_gamma(self,scheme_name='PBE-rel-gamma'):
        """
        Generates calculation scheme for ionic relaxation with PBE functional. Steps: \n
            '1-PBE-SCF-Gamma': Electronic SCF in Gamma point\n
            '2-PBE-OPT-Gamma': Ionic relaxation in Gamma point with maximum 100 ionic steps (NSW=100), Energy change between electronic steps of 1e-05 (EDIFF=1e-05), \n
                         Force convergence criterion of 0.05 eV/Amstrong (EDIFFG=-0.05).
            '3-PBE-SCF': Electronic SCF \n
            '4-PBE-OPT': Ionic relaxation with maximum 100 ionic steps (NSW=100), Energy change between electronic steps of 1e-05 (EDIFF=1e-05), \n
                         Force convergence criterion of 0.05 eV/Amstrong (EDIFFG=-0.05).

        Returns
        -------
        List of VaspJob objects
        """        
        stepnames = ['1-PBE-SCF-Gamma','2-PBE-OPT-Gamma','3-PBE-SCF','4-PBE-OPT']
        jobs = []

        sn = 1
        vaspjob = self.pbe_scf_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 2 
        vaspjob = self.pbe_ionic_rel_gamma(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
              
        sn = 3
        vaspjob = self.pbe_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4 
        vaspjob = self.pbe_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        return jobs
    
    
    def pbe_vol_relaxation(self,scheme_name='PBE-Vrel'):
        """
        Generates calculation scheme for ionic and cell relaxation with PBE functional. Steps: \n
            '1-PBE-SCF': Electronic SC \n
            '2-PBE-OPT': Ionic relaxation with maximum 100 ionic steps (NSW=100), Energy change between electronic steps of 1e-05 (EDIFF=1e-05), \n
                         Force convergence criterion of 0.05 eV/Amstrong (EDIFFG=-0.05) \n
            '3-PBE-VOPT': Same parameters for ionic relaxation and relaxation of the cell (ISIF=3).

        Returns
        -------
        List of VaspJob objects
        """
        stepnames = ['1-PBE-SCF','2-PBE-OPT','3-PBE-VOPT']
        jobs = []
        sn = 1
        vaspjob = self.pbe_scf(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 2 
        vaspjob = self.pbe_ionic_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
 
        sn = 3 
        vaspjob = self.pbe_ionic_vol_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)       
        
        return jobs
                 



        
    