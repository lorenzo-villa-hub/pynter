#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 13:43:48 2025

@author: villa
"""
import numpy as np

from pymatgen.io.vasp.inputs import Kpoints

from pynter.vasp.schemes.core import InputSets


class VariationsSchemes(InputSets):


    def charge_states(self,charges,locpot=True):
        """
        Generate calculation schemes for calculations of different charge states.
        Only the NELECT parameter in INCAR is varied.
        A list of desired charges (float or int) is required.
        
        Returns:
            List of VaspJob objects        
        """  
        jobs = []        
        val = {}
        for p in self.potcar:
            val[p.element] = p.nelectrons        
        nelect = sum([ val[el.symbol]*coeff for el,coeff in self.structure.composition.items()])
        
        for q in charges:
            if isinstance(q,float) and not q.is_integer():
                q = np.around(q,decimals=1)
            stepname = f'q{q}'
            vaspjob = self.get_vaspjob(add_to_job_name=stepname,add_to_job_path=stepname)
            vaspjob.incar['NELECT'] = nelect - q 
            if locpot:
                vaspjob.incar['LVTOT'] = '.TRUE.' #most likely needed for corrections in defect calculations
            jobs.append(vaspjob)
            
        return jobs    


    def convergence_encut(self,encuts=[]):
        """
        Scheme for energy cutoff convergence

        Parameters
        ----------
        encuts : (list), optional
            List of energy cutoffs to use. If not provided a range from 200 to 700 eV is used.The default is [].

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        jobs = []        
        if not encuts:
            encuts = range(300,800,100)

        for ec in encuts:           
            print(ec)
            stepname = f'cutoff{ec}'
            vaspjob = self.pbe_scf(add_to_job_name=stepname,add_to_job_path=stepname)            
            vaspjob.incar['ENCUT'] = ec
            jobs.append(vaspjob)
            
        return jobs
        

    def convergence_kpoints(self,kpoints_meshes=[]):
        """
        Scheme for kpoints convergence.

        Parameters
        ----------
        kpoints_meshes : (list), optional
            List of kpoints meshes to use. If not provided a range from 2x2x2 to 7x7x7 is used.The default is [].

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        jobs = []        
        if not kpoints_meshes:
            kpoints_meshes = []
            for k in range(2,9):
                kpoints_meshes.append((k,k,k))

        for kmesh in kpoints_meshes:
            stepname = 'k%ix%ix%i' %(kmesh[0],kmesh[1],kmesh[2])
            vaspjob = self.get_vaspjob(add_to_job_name=stepname,add_to_job_path=stepname)
            vaspjob.incar.pop('KPAR', None)
            vaspjob.incar['NSW'] = 0
            vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=kmesh)
            jobs.append(vaspjob)
            
        return jobs
            


    def fractional_charge_linearity(self,nelect=None):
        """
        Generate calculation scheme for occupation linearity test.
        The number of electrons are scanned from NELECT to NELECT + 1 with interval of 0.2.
        
        Parameters
        --------
        nelect: (int)
            NELECT to start from. If None the total valence electrons from Potcar are used.
        Returns:
            List of VaspJob objects
        """        
        jobs = []        
        val = {}
        for p in self.potcar:
            val[p.element] = p.nelectrons 
        if not nelect:
            nelect = sum([ val[el]*self.structure.composition.as_dict()[el] for el in self.structure.composition.as_dict()])
        
        for q in np.arange(0,1.2,0.2):
            q = np.around(q,decimals=1)
            stepname = f'q{q}'
            vaspjob = self.get_vaspjob(add_to_job_name=stepname,add_to_job_path=stepname)
            vaspjob.incar['NELECT'] = nelect + q
            jobs.append(vaspjob)
            
        return jobs