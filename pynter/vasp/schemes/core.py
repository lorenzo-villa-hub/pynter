#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 13:30:52 2025

@author: villa
"""

import os
import os.path as op

from pymatgen.io.vasp.inputs import Kpoints, Poscar, Potcar, VaspInput, Incar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

from pynter.hpc.slurm import JobSettings
from pynter.jobs.vasp.vasp_jobs import VaspJob


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
            

class DefaultJobSettings:
    
    def __init__(self,
                 vasp_exe=None,
                 vasp_sbatch=None,
                 modules=None,
                 lines_before_srun=None,
                 lines_after_srun=None):
        """
        Class to handle VASP specific part of sbatch submission script.
        Args that are not provided are taken from vasp.yml file in ~/.pynter.

        Parameters
        ----------
        vasp_exe : (str)
            Path to vasp executable. 
        vasp_sbatch : (dict)
            Sbatch dictionary tailored for vasp.
        modules : (list)
            Name of modules to load with "module load".
        lines_before_srun : (list)
            Script lines to be placed before the "srun" command.
        lines_after_srun : (list)
            Script lines to be placed before the "srun" command.
        """
        
        default_settings = SETTINGS['vasp']['job_settings_default']
        self.vasp_exe = vasp_exe or default_settings['vasp_exe']
        self.vasp_sbatch = vasp_sbatch or default_settings['vasp_sbatch']
        self.modules = modules or default_settings['modules']
        self.lines_before_srun = lines_before_srun or default_settings['lines_before_srun']
        self.lines_after_srun = lines_after_srun or default_settings['lines_after_srun']
    
    @property
    def vasp_script_lines(self):
        if self.modules:
            lines = self.get_modules_lines()
        else:
            lines = []
        if self.lines_before_srun:
            lines += self.lines_before_srun
        if self.vasp_exe:
            lines += [self.get_srun_line()]
        if self.lines_after_srun:
            lines += self.lines_after_srun
        return lines
        
    def get_modules_lines(self):
        if self.modules:
            return [f'module load {module}' for module in self.modules] 

    def get_srun_line(self):
        return 'srun %s' %self.vasp_exe

    def get_updated_job_settings(self,job_settings):
        """
        Get a copy of a JobSettings object updated with VASP specific settings
        """
        job_settings = job_settings.copy()
        for k,v in self.vasp_sbatch.items():
            job_settings[k] = v
                                
        lines = self.vasp_script_lines
        for line_to_add in lines:
            if all([line_to_add not in line for line in job_settings['script_lines']]):
                job_settings['script_lines'].append(line_to_add)

        return job_settings
    
    

class InputSets:
    
    def __init__(self,
                 path=None,
                 vaspinput=None,
                 structure=None,
                 incar_settings=None,
                 kpoints=None,
                 potcar=None,
                 job_settings=None,
                 name=None,
                 add_parent_folder=False,
                 **kwargs):
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
            Dictionary with job settings to create job script, parameters are defined in JobSettings class function. The default is None.\n
            If job_settings is None, the 'job-name' key will be added, the value is the 'name' argument if provided, if 'name' arg is \n
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
            if job_settings:
                self.job_settings = JobSettings(**job_settings)
            elif name:
                self.job_settings = JobSettings(sbatch={'job-name':name})
            else:
                self.job_settings = JobSettings(sbatch={'job-name':'no_name'})
            self.name = name if name != None else self.job_settings['job-name'] 
            
            if self.name: # this return false if name is '', the previuos line considers only if name is None
                self.job_settings['job-name'] = self.name
            
         #   if all(['srun' not in line for line in self.job_settings['script_lines']]):
            self.job_settings = DefaultJobSettings(**kwargs).get_updated_job_settings(self.job_settings)    
            
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


    def get_vaspjob(self,add_to_job_name=None,add_to_job_path=None):
        """
        Generate VaspJob object from the input settings of the class.

        Parameters
        ----------
        add_to_job_name : (str), optional
            String to be added to 'name' key in job_settings dictionary and to name attribute of VaspJob.
        add_to_job_path : (str), optional
            String to be added to self.path. The complete path will be input in 'path' arg of VaspJob.

        Returns
        -------
        VaspJob object
        """        
        incar_settings = self.incar_settings.copy()
        job_settings = self.job_settings.copy()
        
        incar = Incar(incar_settings)
        kpoints = self.kpoints
        poscar = Poscar(self.structure)
        potcar = self.potcar
        vaspinput = VaspInput(incar,kpoints,poscar,potcar)
        
        if add_to_job_name:
            job_settings['job-name'] = '_'.join([job_settings['job-name'],add_to_job_name])        
            jobname = '_'.join([self.name,add_to_job_name])
        else:
            jobname = self.name
            
        if add_to_job_path:
            jobpath = op.join(self.path,add_to_job_path)
        else:
            jobpath = self.path
            
        job_script_filename = job_settings['filename'] if 'filename' in job_settings.keys() else None
        vaspjob = VaspJob(path=jobpath,inputs=vaspinput,job_settings=job_settings,job_script_filename=job_script_filename,name=jobname)
        
        return vaspjob  


    def hse_bs(self,kpoints_bs=None,add_to_job_name='HSE-BS',add_to_job_path='HSE-BS',**kwargs):
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
        vaspjob = self.hse_scf(add_to_job_name,add_to_job_path)
        if kpoints_bs:
            kpoints = kpoints_bs
        else:
            kpoints = DefaultInputs(self.structure).get_kpoints_bs_default(hybrid_mode=True,**kwargs)
        vaspjob.inputs['KPOINTS'] = kpoints
        return vaspjob

    def hse_dos(self,kmesh=3,add_to_job_name='HSE-DOS',add_to_job_path='HSE-DOS'):
        """
        Set up DOS calculation with HSE

        Parameters
        ----------
        kmesh : (Int), optional
            Multiplier for coefficients of the k-mesh for the DOS with respect to first step SCF calculation. The default is 3.
        """
        vaspjob = self.hse_scf(add_to_job_name,add_to_job_path)
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

    def hse_ionic_rel(self,add_to_job_name='HSE-rel',add_to_job_path='HSE-rel'):
        """
        Set up ionic relaxation with HSE
        """
        vaspjob = self.pbe_ionic_rel(add_to_job_name,add_to_job_path)
        vaspjob.incar['LHFCALC'] = '.TRUE.'
        vaspjob.incar['ISYM'] = 3
        vaspjob.job_settings['array_size'] = 7
        return vaspjob

    def hse_ionic_vol_rel(self,add_to_job_name='HSE-rel',add_to_job_path='HSE-rel'):
        """
        Relaxation of atomic positions and cell volume with HSE
        """
        vaspjob = self.hse_ionic_rel(add_to_job_name,add_to_job_path)
        vaspjob.incar['ISIF'] = 3        
        return vaspjob

    def hse_ionic_rel_gamma(self,add_to_job_name='HSE-rel-Gamma',add_to_job_path='HSE-rel-gamma'):
        """
        Set up HSE ionic relaxation only in gamma point 
        """
        vaspjob = self.hse_ionic_rel(add_to_job_name,add_to_job_path)
        vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=(1,1,1))         
        return vaspjob 

    def hse_scf(self,add_to_job_name='HSE-SCF',add_to_job_path='HSE-SCF'):
        """
        Set up SCF calculation with HSE
        """
        vaspjob = self.pbe_scf(add_to_job_name,add_to_job_path)
        vaspjob.incar['LHFCALC'] = '.TRUE.'
        vaspjob.incar['ISYM'] = 3  
        return vaspjob

    def hse_scf_gamma(self,add_to_job_name='HSE-SCF-Gamma',add_to_job_path='HSE-SCF-Gamma'):
        """
        Set up HSE-SCF calculation only in gamma point 
        """
        vaspjob = self.hse_scf(add_to_job_name,add_to_job_path)
        vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=(1,1,1))        
        return vaspjob 

    def hse_standard(self,add_to_job_name='HSE',add_to_job_path='HSE'):
        """
        Set for standard HSE caculation
        """
        vaspjob = self.get_vaspjob(add_to_job_name,add_to_job_path)
        vaspjob.incar['LHFCALC'] = '.TRUE.'
        vaspjob.incar['ISYM'] = 3 
        return vaspjob


    def hubbard(self,ldauu_dict,add_to_job_name,add_to_job_path):
        """
        Set for standard PBE+U calculation. The dict {Element:U_value} is used to construct the INCAR
        """
        vaspjob = self.get_vaspjob(add_to_job_name,add_to_job_path)
        vaspjob.incar['LDAU'] = '.TRUE.'
        vaspjob.incar['LDAUTYPE'] = 2
        vaspjob.incar['LDAUPRINT'] = 2
        vaspjob.incar['LDAUU'] = ' '.join([str(ldauu_dict[el]) for el in ldauu_dict])
        return vaspjob


    def pbe_bs(self,kpoints_bs=None,add_to_job_name='PBE-BS',add_to_job_path='PBE-BS',**kwargs):
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
        vaspjob = self.pbe_scf(add_to_job_name,add_to_job_path)
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

    def pbe_dos(self,kmesh=3,add_to_job_name='PBE-DOS',add_to_job_path='PBE-DOS'):
        """
        Set up DOS calculation with PBE

        Parameters
        ----------
        kmesh : (Int), optional
            Multiplier for coefficients of the k-mesh for the DOS with respect to first step SCF calculation. The default is 3.
        """
        vaspjob = self.pbe_scf(add_to_job_name,add_to_job_path)
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
 
    def pbe_ionic_rel(self,add_to_job_name='PBE-rel',add_to_job_path='PBE-rel'):
        """
        Relaxation of atomic positions with PBE
        """
        vaspjob = self.pbe_scf(add_to_job_name,add_to_job_path)
        vaspjob.incar['NSW'] = 100
        vaspjob.incar['EDIFF'] = 1e-05
        vaspjob.incar['ISIF'] = 2        
        return vaspjob 

    def pbe_ionic_rel_gamma(self,add_to_job_name='PBE-rel-Gamma',add_to_job_path='PBE-rel-gamma'):
        """
        Set up PBE ionic relaxation only in gamma point 
        """
        vaspjob = self.pbe_ionic_rel(add_to_job_name,add_to_job_path)
        vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=(1,1,1))        
        return vaspjob 

    def pbe_ionic_vol_rel(self,add_to_job_name='PBE-rel-vol',add_to_job_path='PBE-rel-vol'):
        """
        Relaxation of atomic positions and cell volume with PBE
        """
        vaspjob = self.pbe_ionic_rel(add_to_job_name,add_to_job_path)
        vaspjob.incar['ISIF'] = 3        
        return vaspjob 
                
    def pbe_scf(self,add_to_job_name='PBE-SCF',add_to_job_path='PBE-SCF'):
        """
        Set up standard input for electronic minimization with PBE
        """        
        vaspjob = self.get_vaspjob(add_to_job_name,add_to_job_path)
        vaspjob.incar['NSW'] = 0
        vaspjob.incar['EDIFF'] = 1e-06        
        return vaspjob      

    def pbe_scf_gamma(self,add_to_job_name='PBE-SCF-Gamma',add_to_job_path='PBE-SCF-Gamma'):
        """
        Set up PBE-SCF calculation only in gamma point 
        """
        vaspjob = self.pbe_scf(add_to_job_name,add_to_job_path)
        vaspjob.inputs['KPOINTS'] = Kpoints().gamma_automatic(kpts=(1,1,1))         
        return vaspjob 