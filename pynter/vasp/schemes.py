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
from pymatgen.io.vasp.inputs import VaspInput, Incar, Poscar, Kpoints

from pynter.vasp.default_inputs import DefaultInputs
from pynter.vasp.jobs import VaspJob, VaspNEBJob
from pynter.defects.defects import create_vacancies, create_substitutions


class InputSets:
    
    def __init__(self,path=None,vaspinput=None,structure=None,incar_settings=None,kpoints=None,potcar=None,
                 job_settings=None,name=None,add_parent_folder=False):
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
        add_parent_folder : (bool), optional
            Add folder to the path names like the name of the InputSets. Default is False.
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
            raise ValueError('You need to provide Structure, either as Poscar object in VaspInput or in "structure" arg')

        if add_parent_folder:
            self.path = op.join(self.path,self.name)

    def __str__(self):
        printout = 'System name:"%s" \n' %self.name
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


    def get_vaspjob(self,setname='',pathname=''):
        """
        Generate VaspJob object from the input settings of the class.

        Parameters
        ----------
        setname : (str), optional
            String to be added to 'name' key in job_settings dictionary and to name attribute of VaspJob. The default is ''.
        pathname : (str), optional
            String to be added to self.path. The complete path will be input in 'path' arg of VaspJob. The default is ''.

        Returns
        -------
        vaspjob : (VaspJob object)
        """        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        job_settings['name'] = '_'.join([self.job_settings['name'],setname])
        
        jobname = '_'.join([self.name,setname])
        jobpath = op.join(self.path,pathname)
        job_script_filename = job_settings['filename'] if 'filename' in job_settings.keys() else None
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,job_script_filename=job_script_filename,name=jobname)
        
        return vaspjob  


    def hse_bs(self,kpoints_bs=None,setname='HSE-BS',pathname='HSE-BS',**kwargs):
        """
        Set up band structure calculation with HSE
        
        Parameters
        ----------
        kpoints_bs : (Pymatgen Kpoints object) , optional
            Pymatgen Kpoints object for the path in BS calculation. The default is None. If None the default high symmetry path \n
            from Pymatgen class HighSymmKpath is obtained from the input Structure with 10 points between high symm k-points.
        kwargs : (dict), optional
            Kwargs for get_kpoints_bs_default in DefaultInputs
        """
        vaspjob = self.hse_scf(setname,pathname)
        if kpoints_bs:
            kpoints = kpoints_bs
        else:
            kpoints = DefaultInputs(self.structure).get_kpoints_bs_default(hybrid_mode=True,**kwargs)
        vaspjob.inputs['KPOINTS'] = kpoints
        return vaspjob

    def hse_dos(self,kmesh=3,setname='HSE-DOS',pathname='HSE-DOS'):
        """
        Set up DOS calculation with HSE

        Parameters
        ----------
        kmesh : (Int), optional
            Multiplier for coefficients of the k-mesh for the DOS with respect to first step SCF calculation. The default is 3.
        """
        vaspjob = self.hse_scf(setname,pathname)
        vaspjob.incar['NEDOS'] = 2000
        vaspjob.incar['ISMEAR'] = -5        
        # multiply by 3 coeff of k-mesh
        kpoints = vaspjob.inputs['KPOINTS']
        style, kpts, kpts_shift = kpoints.style, kpoints.kpts, kpoints.kpts_shift
        kpts_new = []
        kpts_new.append([k*kmesh for k in kpts[0]])
        kpoints = Kpoints(style=style,kpts=kpts_new,kpts_shift=kpts_shift)
        vaspjob.inputs['KPOINTS'] = kpoints
        return vaspjob        

    def hse_ionic_rel(self,setname='HSE-rel',pathname='HSE-rel'):
        """
        Set up ionic relaxation with HSE
        """
        vaspjob = self.pbe_ionic_rel(setname,pathname)
        vaspjob.incar['LHFCALC'] = '.TRUE.'
        vaspjob.incar['ISYM'] = 3
        vaspjob.job_settings['array_size'] = 7
        return vaspjob

    def hse_ionic_vol_rel(self,setname='HSE-rel',pathname='HSE-rel'):
        """
        Relaxation of atomic positions and cell volume with HSE
        """
        vaspjob = self.hse_ionic_rel(setname,pathname)
        vaspjob.incar['ISIF'] = 3        
        return vaspjob

    def hse_ionic_rel_gamma(self,setname='HSE-rel-Gamma',pathname='HSE-rel-gamma'):
        """
        Set up HSE ionic relaxation only in gamma point 
        """
        vaspjob = self.hse_ionic_rel(setname,pathname)
        vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=(1,1,1))         
        return vaspjob 

    def hse_scf(self,setname='HSE-SCF',pathname='HSE-SCF'):
        """
        Set up SCF calculation with HSE
        """
        vaspjob = self.pbe_scf(setname,pathname)
        vaspjob.incar['LHFCALC'] = '.TRUE.'
        vaspjob.incar['ISYM'] = 3  
        vaspjob.job_settings['timelimit'] = '72:00:00'
        vaspjob.job_settings['partition'] = 'test7d'
        return vaspjob

    def hse_scf_gamma(self,setname='HSE-SCF-Gamma',pathname='HSE-SCF-Gamma'):
        """
        Set up HSE-SCF calculation only in gamma point 
        """
        vaspjob = self.hse_scf(setname,pathname)
        vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=(1,1,1))        
        return vaspjob 

    def hse_standard(self,setname='HSE',pathname='HSE'):
        """
        Set for standard HSE caculation
        """
        vaspjob = self.get_vaspjob(setname,pathname)
        vaspjob.incar['LHFCALC'] = '.TRUE.'
        vaspjob.incar['ISYM'] = 3 
        return vaspjob


    def hubbard(self,ldauu_dict,setname,pathname):
        """
        Set for standard PBE+U calculation. The dict {Element:U_value} is used to construct the INCAR
        """
        vaspjob = self.get_vaspjob(setname,pathname)
        vaspjob.incar['LDAU'] = '.TRUE.'
        vaspjob.incar['LDAUTYPE'] = 2
        vaspjob.incar['LDAUPRINT'] = 2
        vaspjob.incar['LDAUU'] = ' '.join([str(ldauu_dict[el]) for el in ldauu_dict])
        return vaspjob


    def pbe_bs(self,kpoints_bs=None,setname='PBE-BS',pathname='PBE-BS',**kwargs):
        """
        Set up band structure calculation with PBE
        
        Parameters
        ----------
        kpoints_bs : (Pymatgen Kpoints object) , optional
            Pymatgen Kpoints object for the path in BS calculation. The default is None. If None the default high symmetry path \n
            from Pymatgen class HighSymmKpath is obtained from the input Structure with 10 points between high symm k-points.
        kwargs : (dict), optional
            Kwargs for get_kpoints_bs_default in DefaultInputs
        """
        vaspjob = self.pbe_scf(setname,pathname)
        vaspjob.incar['ICHARG'] = 11
        vaspjob.incar['LORBIT'] = 11
        vaspjob.incar['ISMEAR'] = 0
        vaspjob.incar['NEDOS'] = 2000
        if kpoints_bs:
            kpoints = kpoints_bs
        else:
            kpoints = DefaultInputs(self.structure).get_kpoints_bs_default(**kwargs)
        vaspjob.inputs['KPOINTS'] = kpoints
        return vaspjob

    def pbe_dos(self,kmesh=3,setname='PBE-DOS',pathname='PBE-DOS'):
        """
        Set up DOS calculation with PBE

        Parameters
        ----------
        kmesh : (Int), optional
            Multiplier for coefficients of the k-mesh for the DOS with respect to first step SCF calculation. The default is 3.
        """
        vaspjob = self.pbe_scf(setname,pathname)
        vaspjob.incar['NEDOS'] = 2000
        vaspjob.incar['ISMEAR'] = -5        
        # multiply by 3 coeff of k-mesh
        kpoints = vaspjob.inputs['KPOINTS']
        style, kpts, kpts_shift = kpoints.style, kpoints.kpts, kpoints.kpts_shift
        kpts_new = []
        kpts_new.append([k*kmesh for k in kpts[0]])
        kpoints = Kpoints(style=style,kpts=kpts_new,kpts_shift=kpts_shift)
        vaspjob.inputs['KPOINTS'] = kpoints
        return vaspjob
 
    def pbe_ionic_rel(self,setname='PBE-rel',pathname='PBE-rel'):
        """
        Relaxation of atomic positions with PBE
        """
        vaspjob = self.pbe_scf(setname,pathname)
        vaspjob.incar['NSW'] = 100
        vaspjob.incar['EDIFF'] = 1e-05
        vaspjob.incar['ISIF'] = 2        
        return vaspjob 

    def pbe_ionic_rel_gamma(self,setname='PBE-rel-Gamma',pathname='PBE-rel-gamma'):
        """
        Set up PBE ionic relaxation only in gamma point 
        """
        vaspjob = self.pbe_ionic_rel(setname,pathname)
        vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=(1,1,1))        
        return vaspjob 

    def pbe_ionic_vol_rel(self,setname='PBE-rel-vol',pathname='PBE-rel-vol'):
        """
        Relaxation of atomic positions and cell volume with PBE
        """
        vaspjob = self.pbe_ionic_rel(setname,pathname)
        vaspjob.incar['ISIF'] = 3        
        return vaspjob 
                
    def pbe_scf(self,setname='PBE-SCF',pathname='PBE-SCF'):
        """
        Set up standard input for electronic minimization with PBE
        """        
        vaspjob = self.get_vaspjob(setname,pathname)
        vaspjob.incar['NSW'] = 0
        vaspjob.incar['EDIFF'] = 1e-06        
        return vaspjob      

    def pbe_scf_gamma(self,setname='PBE-SCF-Gamma',pathname='PBE-SCF-Gamma'):
        """
        Set up PBE-SCF calculation only in gamma point 
        """
        vaspjob = self.pbe_scf(setname,pathname)
        vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=(1,1,1))         
        return vaspjob 

  
class Schemes(InputSets):
    """
    Class to generate and write input files for different calculation schemes in VASP
    """
        
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
            vaspjob = self.get_vaspjob(setname=stepname,pathname=stepname)
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
            stepname = f'cutoff{ec}'
            vaspjob = self.pbe_scf(setname=stepname,pathname=stepname)            
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
            vaspjob = self.get_vaspjob(setname=stepname,pathname=stepname)
            vaspjob.incar.pop('KPAR', None)
            vaspjob.incar['NSW'] = 0
            vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=kmesh)
            jobs.append(vaspjob)
            
        return jobs
            

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
        vaspjob = self.get_vaspjob(setname='eps-el',pathname=scheme_name)    
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
        vaspjob = self.get_vaspjob(setname='eps-ion-lcal',pathname=scheme_name)
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
        vaspjob = self.get_vaspjob(setname='eps-ion-leps',pathname=scheme_name)
        vaspjob.incar['NEDOS'] = 2000
        vaspjob.incar['LEPSILON'] = '.TRUE.'              
        vaspjob.incar['IBRION'] = 8
        vaspjob.incar['NSW'] = 100
 
        return vaspjob


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
            vaspjob = self.get_vaspjob(setname=stepname,pathname=stepname)
            vaspjob.incar['NELECT'] = nelect + q
            jobs.append(vaspjob)
            
        return jobs
        
    
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
        vaspjob = self.pbe_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 2
        vaspjob = self.pbe_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 3 
        vaspjob = self.hse_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4
        vaspjob = self.hse_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
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
        vaspjob = self.pbe_scf_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 2 
        vaspjob = self.pbe_ionic_rel_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 3 
        vaspjob = self.hse_scf_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4 
        vaspjob = self.hse_ionic_rel_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)

        sn = 5 
        vaspjob = self.pbe_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 6 
        vaspjob = self.hse_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 7
        vaspjob = self.hse_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
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
        vaspjob = self.pbe_scf_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 2
        vaspjob = self.pbe_ionic_rel_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 3
        vaspjob = self.hse_scf_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4
        vaspjob = self.hse_ionic_rel_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)

        sn = 5 
        vaspjob = self.pbe_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 6 
        vaspjob = self.pbe_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 7 
        vaspjob = self.hse_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 8 
        vaspjob = self.hse_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
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
        vaspjob = self.hse_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 2
        vaspjob = self.hse_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
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
        vaspjob = self.pbe_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 2
        vaspjob = self.pbe_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        vaspjob.incar['LHFCALC'] = '.FALSE.'
        vaspjob.incar['ISYM'] = 2
        jobs.append(vaspjob)
        
        sn = 3 
        vaspjob = self.hse_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4 
        vaspjob = self.hse_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)

        sn = 5
        vaspjob = self.hse_ionic_vol_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
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
            vaspjob = self.hubbard(ldauu_dict,setname=stepname,pathname=stepname)
            jobs.append(vaspjob)
            
        return jobs

    
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
        vaspjob = self.pbe_ionic_vol_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)  
        
        sn = 2 
        vaspjob = self.pbe_dos(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)

        sn = 3
        vaspjob = self.pbe_bs(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
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
        vaspjob = self.pbe_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 2 
        vaspjob = self.pbe_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
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
        vaspjob = self.pbe_scf_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 2 
        vaspjob = self.pbe_ionic_rel_gamma(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
              
        sn = 3
        vaspjob = self.pbe_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 4 
        vaspjob = self.pbe_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
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
        vaspjob = self.pbe_scf(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
        
        sn = 2 
        vaspjob = self.pbe_ionic_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)
 
        sn = 3 
        vaspjob = self.pbe_ionic_vol_rel(setname=scheme_name+'_%i'%sn ,pathname=stepnames[sn-1])
        jobs.append(vaspjob)       
        
        return jobs
                 

class AdvancedSchemes(Schemes):
    
    
    def defects_pbe_relaxation(self,defects_with_charges,automation=False,locpot=False,rel_scheme='default'):
        """
        Generate jobs for default defect calculation scheme with PBE.

        Parameters
        ----------
        defects_with_charges : (dict)
            List of tuples with Defect object and relative charge state list, for example [(Defect , [q1,q2,q3,...])].
        automation : (bool)
            Add default automation string: 'automation_vasp.py --contcar --chgcar --wavecar --error-check --check-kpoints'.
        locpot : (bool), optional
            Add 'LVTOT=True' in INCAR. The default is False.
        rel_scheme : (str), optional
            Relaxation scheme to use with PBE. 'default for pbe_relaxation and 
            'gamma' pbe_relaxation_gamma. The default is 'default'.

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        jobs = []
        if automation:
            self.job_settings['add_automation'] = 'automation_vasp.py --contcar --chgcar --wavecar --error-check --check-kpoints'

        for df, charge_states in defects_with_charges:
            path = op.join(self.path,df.name)
            defect_name = df.name.replace('(','_').replace(')', '') #replace label parenthesis
            schemes_q = Schemes(path=path,structure=df.defect_structure,incar_settings=self.incar_settings,
                                job_settings=self.job_settings,name=self.name+'_'+defect_name,add_parent_folder=False)
            charge_jobs = schemes_q.charge_states(charge_states,locpot)
            for jq in charge_jobs:
                schemes_rel = Schemes(path=jq.path,vaspinput=jq.inputs,job_settings=jq.job_settings,name=jq.name)
                if rel_scheme == 'default':
                    rel_jobs = schemes_rel.pbe_relaxation()
                elif rel_scheme == 'gamma':
                    rel_jobs = schemes_rel.pbe_relaxation_gamma()
                for jrel in rel_jobs:
                        jobs.append(jrel)
        return jobs


    def interstitials_pbe_relaxation(self):
        """
        Since the interstitials structures neeed to be checked and selected, calculations should be set up manually 
        with defects_pbe_relaxation. If the function needs to be implemented, set up a criterion for the choice of 
        interstitials sites.
        """
        pass

            
    def substitutions_pbe_relaxation(self,elements_to_replace_with_charges,supercell_size=None,automation=False,locpot=False,rel_scheme='default'):
        """
        Generate jobs for default substitutions calculation scheme with PBE.

        Parameters
        ----------
        elements_to_replace_with_charges : (dict)
            Dictionary which indicates the substitution types and the relative charge states list.
            Needs to be in this format:
                {'<new_el>-on-<old_el>':[q0,q1,q2]}, where <new_el> and <old_el> are the symbols of the element to be added and
                the element to be substituted with, respectively. [q0,q1,q2] is the list of corrisponding charge states for this substitution.
            The format of the substitution string is the same that is generated from the function "create_substitutiion_structures".
        supercell_size : (int or numpy array), optional
            Input for the make_supercell function of the Structure class. The default is 1 (the structure is not modified).
        automation : (bool)
            Add default automation string: 'automation_vasp.py --contcar --chgcar --wavecar --error-check --check-kpoints'.
        locpot : (bool), optional
            Add 'LVTOT=True' in INCAR. The default is False.
        rel_scheme : (str), optional
            Relaxation scheme to use with PBE. 'default for pbe_relaxation and 
            'gamma' pbe_relaxation_gamma. The default is 'default'.

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        defects_with_charges = []
        for key,charge_states in elements_to_replace_with_charges.items():
            elements_to_replace = {key.split('-on-')[1] : key.split('-on-')[0]}
            substitution = create_substitutions(self.structure,elements_to_replace=elements_to_replace,
                                                supercell_size=supercell_size)[0]
            defects_with_charges.append((substitution,charge_states))
        jobs = self.defects_pbe_relaxation(defects_with_charges,automation=automation,locpot=locpot,rel_scheme=rel_scheme)
        return jobs        
        

    def vacancies_pbe_relaxation(self,elements_with_charges,supercell_size=None,automation=False,locpot=False,rel_scheme='default'):
        """
        Generate jobs for default vacancies calculation scheme with PBE.

        Parameters
        ----------
        elements_with_charges : (dict)
            Dictionary with element symbols as keys and relative charge states (int) lists as values ({'el':[-1,0,1]}).
        supercell_size : (int or numpy array), optional
            Input for the make_supercell function of the Structure class.
            If None the input structure is not modified. The default is None.
        automation : (bool)
            Add default automation string: 'automation_vasp.py --contcar --chgcar --wavecar --error-check --check-kpoints'.
        locpot : (bool), optional
            Add 'LVTOT=True' in INCAR. The default is False.
        rel_scheme : (str), optional
            Relaxation scheme to use with PBE. 'default for pbe_relaxation and 
            'gamma' pbe_relaxation_gamma. The default is 'default'.

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        defects_with_charges = []
        for el,charge_states in elements_with_charges.items():
            vacancy = create_vacancies(structure=self.structure,elements=[el],supercell_size=supercell_size)[0]
            defects_with_charges.append((vacancy,charge_states))
        jobs = self.defects_pbe_relaxation(defects_with_charges,automation=automation,locpot=locpot,rel_scheme=rel_scheme)
        return jobs

        
    
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
            None the value will be: 'no_name'.
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
        self.job_settings = job_settings if job_settings else ({'name':name} if name else {'name':'NEB'})
        self.name = name if name else None
        
        if 'name' not in self.job_settings.keys():
            self.job_settings['name'] = self.name
        if 'name' in self.job_settings.keys() and self.name:
            self.job_settings['name'] = self.name

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
        scheme_name = scheme_name if scheme_name != None else 'CINEB'
        
        inputs = {}
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()

        incar_settings.update({
            'ISYM': 0,
            'EDIFF': 1e-05,
            'EDIFFG': -0.05,
            'ISIF': 2,
            'NSW' : nionic_steps,
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

          
    def neb_pbe(self,nionic_steps=50,scheme_name=None,stepnames=['NEB']):
        """
        Standard NEB job with PBE. The force convergence is set to 0.1 eV/A, relaxation in done with damped 
        dynamics (IBRION=3), symmetry is turned off (ISYM=0). The maximum number of ionic steps is set to 200.
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
            'NSW' : nionic_steps,
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
                    job_settings['add_automation'] = '(cd ../ && automation_vasp_NEB.py)'
                else:
                    job_settings['add_automation'] = 'automation_vasp.py --chgcar --wavecar'
            job_settings['name'] = '_'.join([self.name,scheme_name,image_name])
            
            jobname = '_'.join([self.name,scheme_name,image_name])
            jobpath = op.join(self.path,stepnames[0],image_name)
            vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,name=jobname)
            jobs.append(vaspjob)
            
        return jobs      
        
        
