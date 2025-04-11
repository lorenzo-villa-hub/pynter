#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 13:38:56 2025

@author: villa
"""

from pynter.vasp.schemes.core import InputSets

  
class ElectronicStructureSchemes(InputSets):
    """
    Class to generate and write input files for different calculation schemes in VASP
    """

    def dielectric_properties_complete(self,scheme_name='dielectric-properties',hybrid=False):
        """
        Set up calculations for electronic and ionic contribution in the dielectric constant.
        If hybrid is set to True the LCALCEPS method is used, else the LEPSILON is used.
        """
        jobs = []
        vaspjob = self.dielectric_properties_electronic()
        jobs.append(vaspjob)
        
        if hybrid:
            vaspjob = self.dielectric_properties_ionic_lcalceps()
        else:
            vaspjob = self.dielectric_properties_ionic_lepsilon()
        jobs.append(vaspjob)
        
        return jobs


    def dielectric_properties_electronic(self,scheme_name='eps-electronic'):
        """
        Set calculation for electronic contribution to the dielectric constant (and also dielectric function).
        Uses 'LOPTICS' method in VASP.
        """
        val = {}
        for p in self.potcar:
            val[p.element] = p.nelectrons
        vaspjob = self.get_vaspjob(add_to_job_name='eps-el',add_to_job_path=scheme_name)    
        vaspjob.incar['NEDOS'] = 2000
        vaspjob.incar['LOPTICS'] = '.TRUE.'              
        nelect = sum([ val[el]*self.structure.composition.as_dict()[el] for el in self.structure.composition.as_dict()]) #number of electrons
        nbands = int((nelect/2 + len(self.structure.sites)/2)*3) # 3*VASP default NBANDS (NELECT/2 +NIONS/2)
        vaspjob.incar['NBANDS'] = nbands
        
        return vaspjob
        

    def dielectric_properties_ionic_lcalceps(self,scheme_name='eps-ionic-lcalceps'):
        """
        Set calculation for ionic contribution to the dielectric constant.
        Uses 'LCALCEPS' method in VASP, combined with 'IBRION=6'. Useful for Hybrid calculations where 'LEPSILON' method does not work.
        """
        vaspjob = self.get_vaspjob(add_to_job_name='eps-ion-lcal',add_to_job_path=scheme_name)
        vaspjob.incar['NEDOS'] = 2000
        vaspjob.incar['LCALCEPS'] = '.TRUE.'              
        vaspjob.incar['IBRION'] = 6
        vaspjob.incar['NSW'] = 100
 
        return vaspjob


    def dielectric_properties_ionic_lepsilon(self,scheme_name='eps-ionic-lepsilon'):
        """
        Set calculation for ionic contribution to the dielectric constant.
        Uses 'LEPSILON' method in VASP, combined with 'IBRION=8'. This method does not work with HSE functionals.
        """
        vaspjob = self.get_vaspjob(add_to_job_name='eps-ion-leps',add_to_job_path=scheme_name)
        vaspjob.incar['NEDOS'] = 2000
        vaspjob.incar['LEPSILON'] = '.TRUE.'              
        vaspjob.incar['IBRION'] = 8
        vaspjob.incar['NSW'] = 100
 
        return vaspjob
    
        
    def pbe_electronic_structure(self,kmesh_dos=3, kpoints_bs=None,scheme_name='PBE-el-str'):
        """
        Generates calculation scheme for electronic structure calculations for PBE (in general non hybrid functionals) Steps: \n
            '1-PBE-relax': Relaxation of atomic positions and cell volume \n
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
        List of VaspJob objects
        """
        stepnames = ['1-PBE-relax','2-PBE-DOS','3-PBE-BS']
        jobs = []
        
        sn = 1 
        vaspjob = self.pbe_ionic_vol_rel(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)  
        
        sn = 2 
        vaspjob = self.pbe_dos(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)

        sn = 3
        vaspjob = self.pbe_bs(add_to_job_name=scheme_name+'_%i'%sn ,add_to_job_path=stepnames[sn-1])
        jobs.append(vaspjob)
        
        return jobs