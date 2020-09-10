#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 14:52:51 2020

@author: villa
"""
import os
import os.path as op
import numpy as np
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import VaspInput, Incar, Poscar, Kpoints, Potcar
from pynter.vasp.default_inputs import DefaultInputs
from pynter.slurm.job_script import ScriptHandler
from pynter.data.jobs import VaspJob, VaspNEBJob
from pynter.data.datasets import Dataset


class Schemes:
    """
    Class to generate and write input files for different calculation schemes in VASP
    """
    
    def __init__(self,path=None,vaspinput=None,structure=None,incar_settings=None,kpoints=None,potcar=None,job_settings=None,name=None):
        """
        Parameters
        ----------
        path : (str)
            Main path for the scheme. If None the current work dir is used. The default is None.
        vaspinput : (Pymatgen VaspInput object), optional
            Set of VASP inputs, the default is None. if provided the other inputs of the class (structure,incar_settings,kpoints,potcar) are not needed 
        structure : (Pymatgen Structure object), optional
            Pymatgen Structure object.
        incar_settings : (Dict), optional
            Dictionary with incar flags. The default is None. If None the default settings for PBE functional from the DefaultInputs class are used.
        kpoints : (Pymatgen Kpoints object), optional
            Pymatgen Kpoints object. The default is None. If None the default settings from the DefaultInputs class are used.
        potcar : (Pymatgen kpoints object), optional
            Pymatgen kpoints object. The default is None. If None the default settings from the DefaultInputs class are used.
        job_settings : (Dict), optional
            Dictionary with job settings to create job script, parameters are defined in ScrpitHandler class function. The default is None.\n
            If job_settings is None, the 'name' key will be added, the value is the 'name' argument if provided, if 'name' arg is \n
            None the value will be: 'no_name'.
        name : (str), optional
            Name for the system to set up scheme for. The default is None.
        """
        self.path = op.abspath(path) if path else os.getcwd()
        if vaspinput:
            structure = vaspinput['POSCAR'].structure
            incar_settings = Incar(vaspinput['INCAR'].copy())
            kpoints = vaspinput['KPOINTS']
            potcar = vaspinput['POTCAR']

        self.structure = structure if structure else None
        if structure:
            self.incar_settings = incar_settings if incar_settings else DefaultInputs(self.structure).get_incar_default()
            self.kpoints = kpoints if kpoints else DefaultInputs(self.structure).get_kpoints_default()
            self.potcar = potcar if potcar else DefaultInputs(self.structure).get_potcar()
            self.job_settings = job_settings if job_settings else ({'name':name} if name else {'name':'no_name'})
            self.name = name if name != None else self.job_settings['name'] 
            
            if self.name: # this return false if name is '', the previuos line considers only if name is None
                self.job_settings['name'] = self.name
                
        else:
            raise ValueError('You need to provide Structure, either as Poscar in VaspInput or in "structure" arg')
                

    def __str__(self):
        printout = 'CalculationScheme object, system name:"%s" \n' %self.name
        printout += 'STRUCTURE:\n'
        printout += self.structure.__str__() + '\n'
        printout += 'INCAR:\n'
        printout += '%s \n' %str(self.incar_settings)
        printout += 'KPOINTS:\n'
        printout += self.kpoints.__str__() + '\n'
        printout += 'POTCAR symbols:\n'
        printout += '%s \n' %(', '.join(self.potcar.symbols))
        printout += 'Job settings:\n'
        printout += '%s' %(str(self.job_settings))
        return printout
    
    def __repr__(self):
        return self.__str__()
        

    def charge_states(self,charges,locpot=True):
        """
        Generate calculation schemes for calculations of different charge states.
        Only the NELECT parameter in INCAR is varied.
        A list of desired charges (float or int) is required.
        
        Returns:
            List of Job objects        
        """
        
        jobs = []
        
        val = {}
        for p in self.potcar:
            val[p.element] = p.nelectrons        
        nelect = sum([ val[el]*self.structure.composition.as_dict()[el] for el in self.structure.composition.as_dict()])
        
        for q in charges:
            
            incar_settings = self.incar_settings.copy()
            job_settings = self.job_settings.copy()
            
            if isinstance(q,float) and not q.is_integer():
                q = np.around(q,decimals=1)
            incar_settings['NELECT'] = nelect - q 
            
            if locpot:
                incar_settings['LVTOT'] = '.TRUE.' #most likely needed for corrections in defect calculations
            incar = Incar(incar_settings)
            kpoints = self.kpoints
            poscar = Poscar(self.structure)
            potcar = self.potcar
            vaspinput = VaspInput(incar,kpoints,poscar,potcar)
            stepname = f'q{q}'
            job_settings['name'] = '_'.join([self.job_settings['name'],stepname])
            
            jobname = '_'.join([self.name,stepname])
            jobpath = op.join(self.path,stepname)
            vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
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
            List of Job objects.
        """
        jobs = []        
        if not encuts:
            encuts = range(300,800,100)

        for ec in encuts:
            
            incar_settings = self.incar_settings.copy()
            job_settings = self.job_settings.copy()
            
            incar_settings['ENCUT'] = ec
            incar_settings['NSW'] = 0
            incar = Incar(incar_settings)
            kpoints = self.kpoints
            poscar = Poscar(self.structure)
            potcar = self.potcar
            vaspinput = VaspInput(incar,kpoints,poscar,potcar)
            stepname = f'cutoff{ec}'
            job_settings['name'] = '_'.join([self.job_settings['name'],stepname])
            
            jobname = '_'.join([self.name,stepname])
            jobpath = op.join(self.path,stepname)
            vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
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
            List of Job objects.
        """
        jobs = []        
        if not kpoints_meshes:
            kpoints_meshes = []
            for k in range(2,9):
                kpoints_meshes.append((k,k,k))

        for kmesh in kpoints_meshes:
            
            incar_settings = self.incar_settings.copy()
            job_settings = self.job_settings.copy()
            
            incar_settings.pop('KPAR', None)
            incar_settings['NSW'] = 0
            incar = Incar(incar_settings)
            kpoints = Kpoints().gamma_automatic(kpts=kmesh)
            poscar = Poscar(self.structure)
            potcar = self.potcar
            vaspinput = VaspInput(incar,kpoints,poscar,potcar)
            stepname = 'k%ix%ix%i' %(kmesh[0],kmesh[1],kmesh[2])
            job_settings['name'] = '_'.join([self.job_settings['name'],stepname])
            
            jobname = '_'.join([self.name,stepname])
            jobpath = op.join(self.path,stepname)
            vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
            jobs.append(vaspjob)
        
        return jobs

               
    def dielectric_properties_electronic(self,scheme_name=None):
        """
        Set calculation for electronic contribution to the dielectric constant (and also dielectric function).
        Uses 'LOPTICS' method in VASP.
        """
        
        scheme_name = scheme_name if scheme_name != None else 'eps-el'
        stepnames = ['eps-electronic']
        jobs = []

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
        
        jobname = '_'.join([self.name,scheme_name,stepnames[0]])
        jobpath = op.join(self.path,stepnames[0])
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        return jobs[0]
        

    def dielectric_properties_ionic_lcalceps(self,scheme_name=None):
        """
        Set calculation for ionic contribution to the dielectric constant.
        Uses 'LCALCEPS' method in VASP, combined with 'IBRION=6'. Useful for Hybrid calculations where 'LEPSILON' method does not work.
        """
        
        scheme_name = scheme_name if scheme_name != None else 'eps-ion-lcal'
        stepnames = ['eps-ionic-lcalceps']
        jobs = []

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
        
        jobname = '_'.join([self.name,scheme_name,stepnames[0]])
        jobpath = op.join(self.path,stepnames[0])
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        return jobs[0]


    def dielectric_properties_ionic_lepsilon(self,scheme_name=None):
        """
        Set calculation for ionic contribution to the dielectric constant.
        Uses 'LEPSILON' method in VASP, combined with 'IBRION=8'. This method does not work with HSE functionals.
        """
        
        scheme_name = scheme_name if scheme_name != None else 'eps-el-leps'
        stepnames = ['eps-ionic-lepsilon']
        jobs = []

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
        
        jobname = '_'.join([self.name,scheme_name,stepnames[0]])
        jobpath = op.join(self.path,stepnames[0])
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        return jobs[0]


    def fractional_charge_linearity(self):
        """
        Generate calculation scheme for occupation linearity test.
        The number of electrons are scanned from NELECT to NELECT + 1 with interval of 0.2.
        
        Returns:
            List of Job objects
        """
        
        jobs = []
        
        val = {}
        for p in self.potcar:
            val[p.element] = p.nelectrons        
        nelect = sum([ val[el]*self.structure.composition.as_dict()[el] for el in self.structure.composition.as_dict()])
        
        for q in np.arange(0,1.2,0.2):
            
            incar_settings = self.incar_settings.copy()
            job_settings = self.job_settings.copy()
            
            q = np.around(q,decimals=1)
            incar_settings['NELECT'] = nelect + q
            
            incar = Incar(incar_settings)
            kpoints = self.kpoints
            poscar = Poscar(self.structure)
            potcar = self.potcar
            vaspinput = VaspInput(incar,kpoints,poscar,potcar)
            stepname = f'q{q}'
            job_settings['name'] = '_'.join([self.job_settings['name'],stepname])
            
            jobname = '_'.join([self.name,stepname])
            jobpath = op.join(self.path,stepname)
            vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
            jobs.append(vaspjob)
            
        return jobs
        
    
    def hse_rel(self,scheme_name=None):
        """
        Generates calculation scheme for ionic relaxation for HSE. Steps: \n
            '1-PBE-SCF': Electronic SCF with PBE \n
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE \n
            '3-HSE-SCF': Electronic SCF with HSE \n
            '4-HSE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE
            
        Returns:
            List of Job objects
        """
        
        scheme_name = scheme_name if scheme_name != None else 'HSE-rel'
        stepnames = ['1-PBE-SCF','2-PBE-OPT','3-HSE-SCF','4-HSE-OPT']
        jobs = []
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        sn = 1 #set step 1        
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 2  # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 3 # set step 3
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 4 # set step 4
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
    
        return jobs
    
    
    def hse_rel_gamma(self,scheme_name=None):
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
            List of Job objects
        """
        
        scheme_name = scheme_name if scheme_name != None else'HSE-rel-gamma'
        stepnames = ['1-PBE-SCF-Gamma','2-PBE-OPT-Gamma','3-HSE-SCF-Gamma',
                  '4-HSE-OPT-Gamma','5-PBE-SCF','6-HSE-SCF','7-HSE-OPT']
        jobs = []
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        sn = 1 #set step 1
        incar_settings['NSW'] = 0
        incar_settings['LHFCALC'] = '.FALSE.'
        incar_settings['ISYM'] = 2
        incar = Incar(incar_settings)
        kpoints = Kpoints().gamma_automatic(kpts=(1,1,1))
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 2 # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 3 # set step 3
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 4 # set step 4
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
    
        sn = 5 # set step 5
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 6 # set step 6
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 7 # set step 7
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        return jobs
    
    
    def hse_rel_gamma_extended(self,scheme_name=None):
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
            List of Job objects
        """

        scheme_name = scheme_name if scheme_name != None else'HSE-rel-gamma-ext'
        stepnames = ['1-PBE-SCF-Gamma','2-PBE-OPT-Gamma','3-HSE-SCF-Gamma',
                  '4-HSE-OPT-Gamma','5-PBE-SCF','6-PBE-OPT','7-HSE-SCF','8-HSE-OPT']
        jobs = []
        
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        sn = 1 #set step 1
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = Kpoints().gamma_automatic(kpts=(1,1,1))
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 2 # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 3 # set step 3
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 4 # set step 4
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
    
        sn = 5 # set step 5
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)

        sn = 6 # set step 6
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 7 # set step 7
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 8 # set step 8
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        return jobs
                      

    def hse_rel_short(self,scheme_name=None):
        """
        Generates calculation scheme for ionic relaxation for HSE (short version i.e. no PBE preliminary calculation). Steps: \n
            '1-PBE-SCF': Electronic SCF with HSE \n
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE

        Returns:
            List of Job objects
        """

        scheme_name = scheme_name if scheme_name != None else 'HSE-rel-short'
        stepnames = ['1-HSE-SCF','2-HSE-OPT']        
        jobs = []
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        sn = 1 # set step 1
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
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 2 # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['nodes'] = 8
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
    
        return jobs 

    
    def hse_vol_rel(self,scheme_name=None):
        """
        Generates calculation scheme for structure relaxations for HSE including cell volume relaxation. Steps: \n
            '1-PBE-SCF': Electronic SCF with PBE
            '2-PBE-OPT': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with PBE
            '3-HSE-SCF': Electronic SCF with HSE
            '4-HSE-OPT-Gamma': Relaxation of atomic positions with ISIF = 2 and EDIFFG = 0.05 eV/A with HSE
            '5-HSE-VOPT': Cell volume relaxation and ionic relaxation with HSE (ISIF=3)

        Returns:
            List of Job objects
        """
        
        scheme_name = scheme_name if scheme_name != None else 'HSE-Vrel'
        stepnames = ['1-PBE-SCF','2-PBE-OPT','3-HSE-SCF','4-HSE-OPT','5-HSE-VOPT']
        jobs = []
        
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        sn = 1 #set step 1
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 2
        incar_settings['LHFCALC'] = '.FALSE.'
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 2 # set step 2
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 3 # set step 3
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-06
        incar_settings['NSW'] = 0
        incar_settings['ISYM'] = 3
        incar_settings['LHFCALC'] = '.TRUE.'
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 4  # set step 4
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['array_size'] = 7
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 5 # set volume relaxation
        job_settings = self.job_settings.copy()
        incar_settings['ISIF'] = 3
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)

        return jobs
    
    
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
        List of Job objects
        """
        
        scheme_name = scheme_name if scheme_name != None else 'U_tuning'
        jobs = []
        
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
            stepname = f'U_{u}'
            job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,stepname])
            
            jobname = '_'.join([self.name,scheme_name,stepname])
            jobpath = op.join(self.path,scheme_name,stepname)
            vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
            jobs.append(vaspjob)
            
        return jobs

    
    def pbe_electronic_structure(self,kmesh_dos=3, kpoints_bs=None,scheme_name=None):
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
        List of Job objects
        """
     
        scheme_name = scheme_name if scheme_name != None else 'PBE-el-str'
        stepnames = ['1-PBE-relax','2-PBE-DOS','3-PBE-BS']
        jobs = []
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        sn = 1 #set step 1 : relaxation
        incar_settings['NSW'] = 100
        incar_settings['ISIF'] = 3
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 2 # set step 2 : DOS calculation
        job_settings = self.job_settings.copy()
        incar_settings['NSW'] = 0
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
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)

        sn = 3 # set step 3 : BS calculation
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
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        return jobs
                 
    
    def pbe_rel(self,scheme_name=None):
        """
        Generates calculation scheme for ionic relaxation with PBE functional. Steps: \n
            '1-PBE-SCF': Electronic SCF \n
            '2-PBE-OPT': Ionic relaxation with maximum 100 ionic steps (NSW=100), Energy change between electronic steps of 1e-05 (EDIFF=1e-05), \n
                         Force convergence criterion of 0.05 eV/Amstrong (EDIFFG=-0.05).

        Returns
        -------
        List of Job objects
        """        

        scheme_name = scheme_name if scheme_name != None else 'PBE-rel'
        stepnames = ['1-PBE-SCF','2-PBE-OPT']
        jobs = []
               
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        sn = 1 #set step 1
        incar_settings['NSW'] = 0
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 2 # set step 2 : ionic relaxation
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        return jobs
    
    
    def pbe_vol_rel(self,scheme_name=None):
        """
        Generates calculation scheme for ionic and cell relaxation with PBE functional. Steps: \n
            '1-PBE-SCF': Electronic SC \n
            '2-PBE-OPT': Ionic relaxation with maximum 100 ionic steps (NSW=100), Energy change between electronic steps of 1e-05 (EDIFF=1e-05), \n
                         Force convergence criterion of 0.05 eV/Amstrong (EDIFFG=-0.05) \n
            '3-PBE-VOPT': Same parameters for ionic relaxation and relaxation of the cell (ISIF=3).

        Returns
        -------
        List of Job objects
        """
        
        scheme_name = scheme_name if scheme_name != None else 'PBE-Vrel'
        stepnames = ['1-PBE-SCF','2-PBE-OPT','3-PBE-VOPT']
        jobs = []
        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        sn = 1 #set step 1
        incar_settings['NSW'] = 0
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 2 # set step 2 :ionic relaxation
        job_settings = self.job_settings.copy()
        incar_settings['EDIFF'] = 1e-05
        incar_settings['NSW'] = 100
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        sn = 3 # set step 3: volume relaxation
        job_settings = self.job_settings.copy()
        incar_settings['ISIF'] = 3
        incar = Incar(incar_settings)
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],scheme_name,str(sn)])
        
        stepname = stepnames[sn-1]
        jobname = '_'.join([self.name,scheme_name,str(sn)])
        jobpath = op.join(self.path,scheme_name,stepname)
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
        jobs.append(vaspjob)
        
        return jobs
                 
        
class NEBSchemes:
    
    def __init__(self,path,structures,incar_settings=None,kpoints=None,potcar=None,job_settings=None,name=None):
        """
        Parameters
        ----------
        path : (str)
            Main path for the scheme.        
        structures : (list of Pymatgen Structure objects), optional
            List of structures for reaction path.
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
        
        self.path = path
        self.structures = structures
        self.incar_settings = incar_settings if incar_settings else DefaultInputs(self.structures[0]).get_incar_default()
        self.kpoints = kpoints if kpoints else DefaultInputs(self.structures[0]).get_kpoints_default()
        self.potcar = potcar if potcar else DefaultInputs(self.structures[0]).get_potcar()
        self.job_settings = job_settings if job_settings else ({'name':name} if name else {'name':'NEB'})
        self.name = name if name else None
        
        if 'name' not in self.job_settings.keys():
            self.job_settings['name'] = self.name
        if 'name' in self.job_settings.keys() and self.name:
            self.job_settings['name'] = self.name


        for key in ['KPAR','NCORE','NPAR']:
            if key in self.incar_settings:
                del self.incar_settings[key]

        self.images = len(self.structures)-2


    def cineb_pbe(self,scheme_name=None,stepnames=['CINEB']):
        """
        Climbing image NEB job with PBE. The force convergence is set to 0.05 eV/A, relaxation in done with damped 
        dynamics (IBRION=3), symmetry is turned off (ISYM=0). The maximum number of ionic steps is set to 25.
        """
        scheme_name = scheme_name if scheme_name != None else 'CINEB'
        
        inputs = {}
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()

        incar_settings.update({
            'ISYM': 0,
            'EDIFF': 1e-05,
            'EDIFFG': -0.05,
            'ISIF': 2,
            'NSW' : 25,
            'IMAGES' : self.images,
            'SPRING' : -5,
            'IOPT' : 0,
            'IBRION' : 3,
            'POTIM' : 0.05,
            'ICHAIN' : 0,
            'LCLIMB' : '.TRUE.',
            'LTANGENTOLD' : '.FALSE.',
            'LDNEB' : '.FALSE.',
            'LNEBCELL' : '.FALSE.'
            })

        incar = Incar(incar_settings)
        kpoints = self.kpoints
        potcar = self.potcar
        
        inputs['structures'] = self.structures
        inputs['INCAR'] = incar
        inputs['KPOINTS'] = kpoints
        inputs['POTCAR'] = potcar
        
        job_settings['nodes'] = self.images
        job_settings['path_exe'] = '/home/lv51dypu/vasp.5.4.4_vtstcode/bin/vasp_std'
        if 'add_automation' not in job_settings:
            job_settings['add_automation'] = None        
        job_settings['name'] = '_'.join([self.name,scheme_name])


        jobname = '_'.join([self.name,scheme_name])
        jobpath = op.join(self.path,stepnames[0])
        nebjob = VaspNEBJob(path=jobpath,inputs=inputs,job_settings=job_settings,name=jobname)

        return nebjob

 
    def neb_complete_pbe(self,scheme_name=None):
        """
        NEB scheme with 3 steps:
            - Step 1: preconvergence (SCF) of the images with EDIFF = 1e-04
            - Step 2: NEB calculation with NSW = 25
            - Step 3: CINEB with NSW = 25
        """
        scheme_name = scheme_name if scheme_name != None else 'NEB_comp'
        jobs = []
        step1 = self.preconverge(scheme_name='NEB-SCF_1',stepnames=['1-SCF'])
        for j in step1:
            jobs.append(j)
        
        step2 = self.neb_pbe(scheme_name='NEB_2',stepnames=['2-NEB'])
        jobs.append(step2)
        
        step3 = self.cineb_pbe(scheme_name='NEB_3',stepnames=['3-CINEB'])
        jobs.append(step3)

        return jobs
 

    def neb_complete_4step_pbe(self,scheme_name=None):
        """
        NEB scheme with 3 steps:
            - Step 1: preconvergence (SCF) of the images with EDIFF = 1e-04
            - Step 2: NEB calculation with NSW = 10
            - Step 3: NEB calculation with NSW = 25
            - Step 4: CINEB with NSW = 25
        """
        scheme_name = scheme_name if scheme_name != None else 'NEB_comp'
        jobs = []
        step1 = self.preconverge(scheme_name='NEB-SCF_1',stepnames=['1-SCF'])
        for j in step1:
            jobs.append(j)
        
        step2 = self.neb_pbe(scheme_name='NEB_2',stepnames=['2-NEB-NSW10'])
        step2.inputs['INCAR']['NSW'] = 10
        jobs.append(step2)

        step3 = self.neb_pbe(scheme_name='NEB_3',stepnames=['3-NEB-NSW25'])
        jobs.append(step3)
        
        step4 = self.cineb_pbe(scheme_name='NEB_4',stepnames=['4-CINEB'])
        jobs.append(step4)

        return jobs

          
    def neb_pbe(self,scheme_name=None,stepnames=['NEB']):
        """
        Standard NEB job with PBE. The force convergence is set to 0.1 eV/A, relaxation in done with damped 
        dynamics (IBRION=3), symmetry is turned off (ISYM=0). The maximum number of ionic steps is set to 25.
        """
        scheme_name = scheme_name if scheme_name != None else 'NEB'
        
        inputs = {}
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()

        incar_settings.update({
            'ISYM': 0,
            'EDIFF': 1e-04,
            'EDIFFG': -0.1,
            'ISIF': 2,
            'NSW' : 25,
            'IMAGES' : self.images,
            'SPRING' : -5,
            'IOPT' : 0,
            'IBRION' : 3,
            'POTIM' : 0.05,
            'ICHAIN' : 0,
            'LCLIMB' : '.FALSE.',
            'LTANGENTOLD' : '.FALSE.',
            'LDNEB' : '.FALSE.',
            'LNEBCELL' : '.FALSE.'
            })

        incar = Incar(incar_settings)
        kpoints = self.kpoints
        potcar = self.potcar
        
        inputs['structures'] = self.structures
        inputs['INCAR'] = incar
        inputs['KPOINTS'] = kpoints
        inputs['POTCAR'] = potcar
        
        job_settings['nodes'] = self.images
        if 'add_automation' not in job_settings:
            job_settings['add_automation'] = 'automation_vasp_NEB.py'        
        job_settings['name'] = '_'.join([self.name,scheme_name])


        jobname = '_'.join([self.name,scheme_name])
        jobpath = op.join(self.path,stepnames[0])
        nebjob = VaspNEBJob(path=jobpath,inputs=inputs,job_settings=job_settings,name=jobname)

        return nebjob


    def preconverge(self,scheme_name=None,stepnames=['NEB-preconverge']):
        """
        Do SCF convergence for all images before starting NEB calculation
        """      
        scheme_name = scheme_name if scheme_name!=None else 'SCF'

        self.incar_settings.update({
            'EDIFF' : 1e-04,
            'NSW': 0
            })

        jobs = []
        structures = self.structures
        for s in structures:
            
            index = structures.index(s)
            image_name = str(index).zfill(2)
            
            incar = Incar(self.incar_settings)
            kpoints = self.kpoints
            potcar = self.potcar
            poscar = Poscar(s)
            vaspinput = VaspInput(incar, kpoints, poscar, potcar)
            
            job_settings = self.job_settings.copy()
            if 'add_automation' not in job_settings:
                if index ==  structures.index(structures[-1]):
                    job_settings['add_automation'] = 'automation_vasp_NEB.py'
                else:
                    job_settings['add_automation'] = 'automation_vasp.py --chgcar --wavecar'
            job_settings['name'] = '_'.join([self.name,scheme_name,image_name])
            
            jobname = '_'.join([self.name,scheme_name,image_name])
            jobpath = op.join(self.path,stepnames[0],image_name)
            vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
            jobs.append(vaspjob)
            
        return jobs      
        
        
