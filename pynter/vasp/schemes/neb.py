#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 13:35:52 2025

@author: villa
"""
import os.path as op

from pymatgen.io.vasp.inputs import Poscar, VaspInput, Incar

from pynter.hpc.slurm import JobSettings
from pynter.jobs.vasp.vasp_jobs import VaspJob, VaspNEBJob
from pynter.vasp.schemes.core import DefaultInputs


class NEBSchemes:
    
    def __init__(self,path,structures,incar_settings=None,kpoints=None,potcar=None,job_settings=None,name=None,add_parent_folder=False):
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
            None the value will be: 'no_name'. If automations are not desired, set job_settings['add_automation'] equal to ''.
        name : (str), optional
            Name for the system to set up scheme for. The default is None.
        add_parent_folder : (bool), optional
            Add folder to the path names like the name of the InputSets. Default is False.
        """
        
        self.path = path
        self.structures = structures
        self.incar_settings = incar_settings if incar_settings else DefaultInputs(self.structures[0]).get_incar_default()
        self.kpoints = kpoints if kpoints else DefaultInputs(self.structures[0]).get_kpoints_default()
        self.potcar = potcar if potcar else DefaultInputs(self.structures[0]).get_potcar()
        if job_settings:
            self.job_settings = JobSettings(**job_settings)
        elif name:
            self.job_settings = JobSettings(slurm={'job-name':name})
        else:
            self.job_settings = JobSettings(slurm={'job-name':'no_name'})
        self.name = name if name != None else self.job_settings['job-name'] 
        
        if self.name: # this return false if name is '', the previuos line considers only if name is None
            self.job_settings['job-name'] = self.name

        if add_parent_folder:
            self.path = op.join(self.path,self.name)

        for key in ['KPAR','NCORE','NPAR']:
            if key in self.incar_settings:
                del self.incar_settings[key]

        self.images = len(self.structures)-2


    def cineb_pbe(self,nionic_steps=50,scheme_name=None,stepnames=['CINEB']):
        """
        Climbing image NEB job with PBE. The force convergence is set to 0.05 eV/A, relaxation is done with damped 
        dynamics (IBRION=3), symmetry is turned off (ISYM=0). The maximum number of ionic steps is set to nionic_steps.
        """
        if not scheme_name:
            scheme_name = stepnames[0]

        nebjob = self.neb_generic_pbe(scheme_name=scheme_name,stepnames=stepnames)
        nebjob.incar.update({
            'ISYM': 0,
            'EDIFF': 1e-05,
            'EDIFFG': -0.05,
            'ISIF': 2,
            'NSW' : nionic_steps,
            'IOPT' : 0,
            'IBRION' : 3,
            'POTIM' : 0.05,
            'LCLIMB':'.TRUE.'
            })

        return nebjob

 
    def neb_complete_pbe(self,nionic_steps=50,scheme_name=None):
        """
        NEB scheme with 3 steps:
            - Step 1: preconvergence (SCF) of the images with EDIFF = 1e-04
            - Step 2: NEB calculation with NSW = nionic_steps
            - Step 3: CINEB with NSW = nionic_steps
        """
        scheme_name = scheme_name if scheme_name != None else 'NEB_comp'
        jobs = []
        step1 = self.preconverge(scheme_name='NEB-SCF_1',stepnames=['1-SCF'])
        for j in step1:
            jobs.append(j)
        
        step2 = self.neb_pbe(nionic_steps=nionic_steps,scheme_name='NEB_2',stepnames=['2-NEB'])
        jobs.append(step2)
        
        step3 = self.cineb_pbe(nionic_steps=nionic_steps,scheme_name='NEB_3',stepnames=['3-CINEB'])
        jobs.append(step3)

        return jobs
 

    def neb_complete_4step_pbe(self,nionic_steps=50,scheme_name=None):
        """
        NEB scheme with 4 steps:
            - Step 1: preconvergence (SCF) of the images with EDIFF = 1e-04
            - Step 2: NEB calculation with NSW = 10
            - Step 3: NEB calculation with NSW = nionic_steps
            - Step 4: CINEB with NSW = nionic_steps
        """
        scheme_name = scheme_name if scheme_name != None else 'NEB_comp'
        jobs = []
        step1 = self.preconverge(scheme_name='NEB-SCF_1',stepnames=['1-SCF'])
        for j in step1:
            jobs.append(j)
        
        step2 = self.neb_pbe(scheme_name='NEB_2',stepnames=['2-NEB-NSW10'])
        step2.inputs['INCAR']['NSW'] = 10
        jobs.append(step2)

        step3 = self.neb_pbe(nionic_steps=nionic_steps,scheme_name='NEB_3',stepnames=['3-NEB'])
        jobs.append(step3)
        
        step4 = self.cineb_pbe(nionic_steps=nionic_steps,scheme_name='NEB_4',stepnames=['4-CINEB'])
        jobs.append(step4)

        return jobs


    def neb_generic_pbe(self,scheme_name=None,stepnames=['NEB-generic']):
        """
        Returns a VaspNEBJon with just necessary changes to INCAR. No particular calculation scheme is set.
        """
        if not scheme_name:
            scheme_name = stepnames[0]
        inputs = {}
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()

        incar_settings.update({
            'IMAGES' : self.images,
            'SPRING' : -5,
            'IOPT' : 0, #default
            'ICHAIN' : 0, #default 
            'LCLIMB' : '.FALSE.', 
            'LTANGENTOLD' : '.FALSE.', #default
            'LDNEB' : '.FALSE.', #default
            'LNEBCELL' : '.FALSE.' #default
            })

        incar = Incar(incar_settings)
        kpoints = self.kpoints
        potcar = self.potcar
        
        inputs['structures'] = self.structures
        inputs['INCAR'] = incar
        inputs['KPOINTS'] = kpoints
        inputs['POTCAR'] = potcar        
        job_settings['job-name'] = '_'.join([self.name,scheme_name])

        jobname = '_'.join([self.name,scheme_name])
        jobpath = op.join(self.path,stepnames[0])
        nebjob = VaspNEBJob(path=jobpath,inputs=inputs,job_settings=job_settings,name=jobname)

        return nebjob       
    

    def neb_pbe(self,nionic_steps=50,scheme_name=None,stepnames=['NEB']):
        """
        Standard NEB job with PBE. The force convergence is set to 0.1 eV/A, relaxation in done with damped 
        dynamics (IBRION=3), symmetry is turned off (ISYM=0). The maximum number of ionic steps is set to 200.
        """
        if not scheme_name:
            scheme_name = stepnames[0]

        nebjob = self.neb_generic_pbe(scheme_name=scheme_name,stepnames=stepnames)
        nebjob.incar.update({
            'ISYM': 0,
            'EDIFF': 1e-04,
            'EDIFFG': -0.1,
            'ISIF': 2,
            'NSW' : nionic_steps,
            'IBRION' : 3,
            'POTIM' : 0.05,
            })

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
            if 'add_automation' not in job_settings  or job_settings['add_automation'] == None:
                if index ==  structures.index(structures[-1]):
                    job_settings['add_automation'] = '(cd ../ && pynter automation vasp-NEB)'
                else:
                    job_settings['add_automation'] = 'pynter automation vasp --chgcar --wavecar'
            job_settings['job-name'] = '_'.join([self.name,scheme_name,image_name])
            
            jobname = '_'.join([self.name,scheme_name,image_name])
            jobpath = op.join(self.path,stepnames[0],image_name)
            vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
            jobs.append(vaspjob)
            
        return jobs      
        