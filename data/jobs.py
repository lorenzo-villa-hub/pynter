
from abc import ABCMeta, abstractmethod
import os
import os.path as op
import shutil
from pymatgen.io.vasp.inputs import VaspInput, Poscar, Incar, Kpoints, Potcar
from pymatgen.io.vasp.outputs import Vasprun, Oszicar
from pymatgen.analysis.transition_state import NEBAnalysis
from pynter.slurm.job_script import ScriptHandler
from pynter.slurm.interface import HPCInterface
from pynter.tools.grep import grep_list
import importlib


class Job:
    
    def __init__(self,path=None,inputs=None,job_settings=None,outputs=None,job_script_filename='job.sh',name=None):
        """
        Class to control and organize inputs and outputs of a generic job.

        Parameters
        ----------
        path : (str), optional
            Path where job is stored. The default is None. If None the work dir is used.
        inputs : (dict), optional
            Dictionary with input data. The default is None.
        job_settings : (dict), optional
            Dictionary with job settings. The default is None. Documentation in ScriptHandler class in slurm.job_script module
        outputs : (dict), optional
            Dictionary with output data. The default is None.
        job_script_filename : (str), optional
            Filename of job script. The default is 'job.sh'.
        name : (str)
            Name of the job. If none the name is searched in the job script.

        """
                
        self.path = path if path else os.getcwd()
        self.inputs = inputs
        self.job_settings = job_settings
        self.outputs = outputs
        self.job_script_filename = job_script_filename
        
        self._localdir = HPCInterface().localdir
        self._workdir = HPCInterface().workdir
        self._path_relative = op.abspath(self.path).replace(self._localdir,'')
        
        self.path_in_hpc = self._workdir + self._path_relative
        
        if name:
            self.name = name
        elif self.job_settings:
            self.name = self.job_settings['name']
        elif op.isfile(op.join(self.path,self.job_script_filename)):
            s = ScriptHandler.from_file(self.path,filename=self.job_script_filename)
            self.name = s.settings['name']
        else:
            self.name = 'no_name'
            
        if not self.job_settings:
            self.job_settings = {}
        self.job_settings['name'] = self.name


    def __str__(self):
        jobclass = self.jobclass
        if hasattr(self,'group'):
            if self.group != '':
                printout = '%s "%s" of group "%s"' %(jobclass, self.name, self.group)
        else:
            self.group = ''
            printout = '%s "%s"' %(jobclass, self.name)
        
        return printout
    
    def __repr__(self):
        return self.__str__()
        
    @property
    def jobclass(self):
        return self.__class__.__name__

        
    def cancel_job(self):
        """Cancel job on HPC"""
        hpc = HPCInterface()
        job_id = self.job_id()
        hpc.cancel_jobs(job_id)
        
        return 

    
    def delete_job_files(self,safety=True):
        """
        Delete Job folder (self.path)

        Parameters
        ----------
        safety : (bool), optional
            Ask confirmation to delete job. The default is True.
        """
        if safety:
            inp = input('Are you sure you want to delete Job %s ?: (y/n)' %self.name)
            if inp in ('y','Y'):
                shutil.rmtree(self.path)
                print('Deleted Job %s'%self.name)
            else:
                print('Job %s is left unchanged'%self.name)
        else:
            shutil.rmtree(self.path)
            print('Deleted Job %s'%self.name)
        return


    @abstractmethod
    def get_inputs(self):
        pass

    @abstractmethod
    def get_outputs(self):
        pass


    def job_id(self):
        """Get job ID from the queue on HPC"""        
        hpc = HPCInterface()
        stdout,stderr = hpc.qstat(printout=False)
        queue = stdout.splitlines()
        job_lines = grep_list(self.name,queue)
        if job_lines == []:
            raise ValueError (f'Job named "{self.name}" is not currently running or pending')
        elif len(job_lines) > 1:
            raise ValueError (f'More than one job named "{self.name}" has been found in queue:\n{stdout}')
        else:
            job_line = job_lines[0].split()
            job_id = job_line[0]
        
        return job_id


    def job_queue(self):
        """
        Print job queue from HPC on screen
        
        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        hpc = HPCInterface()
        stdout,stderr = hpc.qstat()
        
        return stdout,stderr
        

    def run_job(self,write_input=True,sync=True):
        """
        Run job on HPC. Input files are automatically written and sync to HPC is performed.

        Parameters
        ----------
        write_input : (bool), optional
            Write input file stored in "inputs" dictionary. The default is True.
        sync : (bool), optional
            Sync files to HPC before running. The default is True

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        if write_input:
            self.write_input()
        hpc = HPCInterface()
        if sync:
            self.sync_to_hpc()
        stdout,stderr = hpc.sbatch(path=self.path_in_hpc,job_script_filename=self.job_script_filename)
        
        return stdout,stderr
    

    def sync_from_hpc(self,stdouts=False):
        """
        Sync job data from HPC to local machine

        Parameters
        ----------
        stdouts : (bool), optional
            Return output and error strings. The default is False.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.

        """
        hpc = HPCInterface()
        abs_path = op.abspath(self.path)
        localdir = abs_path 
        stdout,stderr = hpc.rsync_from_hpc(remotedir=self.path_in_hpc,localdir=localdir)
        if stdouts:
            return stdout,stderr
        else:
            return
        
        
    def sync_to_hpc(self,stdouts=False):
        """
        Sync job data from local machine to HPC

        Parameters
        ----------
        stdouts : (bool), optional
            Return output and error strings. The default is False.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.

        """
        hpc = HPCInterface()
        abs_path = op.abspath(self.path)
        localdir = abs_path 
        stdout,stderr = hpc.rsync_to_hpc(localdir=localdir,remotedir=self.path_in_hpc)
        if stdouts:
            return stdout,stderr
        else:
            return


    def status(self):
        """
        Get job status from HPC. 

        Returns
        -------
        status : (str)
            Job status. Possible status are 'PENDING','RUNNING','NOT IN QUEUE'.
        """
        hpc = HPCInterface()
        stdout,stderr = hpc.qstat(printout=False)
        queue = stdout.splitlines()
        job_lines = grep_list(self.name,queue)
        if job_lines == []:
            status = 'NOT IN QUEUE'
        elif len(job_lines) > 1:
            raise ValueError (f'More than one job named "{self.name}" has been found in queue:\n{stdout}')
        else:
            job_line = job_lines[0].split()
            status = job_line[4]
            if status == 'PD':
                status = 'PENDING'
            if status == 'R':
                status = 'RUNNING'
            if status == 'CG':
                status = 'COMPLETED'
            
        return status
 
           
    @abstractmethod
    def write_input():
        pass
        
     
class VaspJob(Job):
 
    
    def as_dict(self):
        """
        Returns:
            Json-serializable dict representation of VaspJob
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "path": self.path,
             "inputs": self.inputs.as_dict(),
             "job_settings": self.job_settings,
             "outputs": {key: value.as_dict() for key,value in self.outputs.items()},
             "job_script_filename":self.job_script_filename,
             "name":self.name}
        return d
        
      
    @staticmethod
    def from_directory(path,job_script_filename='job.sh'):
        """
        Builds VaspJob object from data stored in a directory. Input files are read using Pymatgen VaspInput class.
        Output files are read usign Pymatgen Vasprun class.
        Job settings are read from the job script file.

        Parameters
        ----------
        path : (str)
            Path were job data is stored.
        job_script_filename : (str), optional
            Filename of job script. The default is 'job.sh'.

        Returns
        -------
        VaspJob object.
        
        """
                
        inputs = VaspInput.from_directory(path)
        outputs = {}
        if op.isfile(op.join(path,'vasprun.xml')):
            try:
                outputs['Vasprun'] = Vasprun(op.join(path,'vasprun.xml'))
            except:
                print('Warning: Reading of vasprun.xml in "%s" failed'%path)
                outputs['Vasprun'] = None
        
        s = ScriptHandler.from_file(path,filename=job_script_filename)
        job_settings =  s.settings
        
        return VaspJob(path,inputs,job_settings,outputs)


    def delete_outputs(self,safety=True):
        """
        Delete files that aren't input files (INCAR,KPOINTS,POSCAR,POTCAR)
        """
        if safety:
            inp = input('Are you sure you want to delete outputs of Job %s ?: (y/n)' %self.name)
            if inp in ('y','Y'):
                delete = True
            else:
                delete = False
        else:
            delete= True
            
        if delete:                
            files = [f for f in os.listdir(self.path) if os.path.isfile(os.path.join(self.path, f))]
            for f in files:
                if f not in ['INCAR','KPOINTS','POSCAR','POTCAR',self.job_script_filename]:
                    os.remove(os.path.join(self.path,f))
                    print('Deleted file %s'%os.path.join(self.path,f))   
        return


    @property
    def incar(self):
        return self.inputs['INCAR']
    
    @property
    def kpoints(self):
        return self.inputs['KPOINTS']
    
    @property
    def poscar(self):
        return self.inputs['POSCAR']
    
    @property
    def potcar(self):
        return self.inputs['POTCAR']

    @property
    def vasprun(self):
        if 'Vasprun' in self.outputs.keys():
            return self.outputs['Vasprun']
        else:
            raise ValueError('"vasprun.xml" file is not present in Job directory')
    
    @property
    def formula(self):
        """Complete formula from initial structure (read with Pymatgen)"""
        if self.initial_structure:
            return self.initial_structure.composition.formula
        else:
            return None

        
    @property
    def initial_structure(self):
        """Initial structure read vasprun.xml if job is converged, otherwise read from Poscar file with Pymatgen"""
        poscar_file = op.join(self.path,'POSCAR')
        if self.is_converged:
            return self.vasprun.structures[0]
        elif op.isfile(poscar_file):
            return Poscar.from_file(poscar_file).structure
        elif self.poscar:
            poscar = self.poscar
            return poscar.structure            
        else:
            return None
    
    
    @property
    def is_converged(self):
        """
        Reads Pymatgen Vasprun object and returns "True" if the calculation is converged,
        "False" if reading failed, and "None" if is not present in the outputs dictionary.
        """
        is_converged = None
        if self.outputs:
            if 'Vasprun' in self.outputs.keys():
                is_converged = False
                if self.vasprun:
                    vasprun = self.vasprun
                    conv_el, conv_ionic = False, False
                    if vasprun:
                        conv_el = vasprun.converged_electronic
                        conv_ionic = vasprun.converged_ionic
                    if conv_el and conv_ionic:
                        is_converged = True
                    
        return is_converged
    
    
    def charge(self):
        """
        Charge of the system calculated as the difference between the value of "NELECT"
        in the INCAR and the number of electrons in POTCAR. If "NELECT" is not present 
        charge is set to 0.
        """
        charge = 0
        if 'NELECT' in self.incar.keys():
            nelect = self.incar['NELECT']
            val = {}
            for p in self.potcar:
                val[p.element] = p.nelectrons
            neutral = sum([ val[el.symbol]*self.initial_structure.composition[el] 
                           for el in self.initial_structure.composition])
            charge = neutral - nelect
        return charge


    def energy_gap(self):
        """Energy gap read from vasprun.xml with Pymatgen"""
        vasprun = self.vasprun
        band_gap = vasprun.eigenvalue_band_properties[0]
        return band_gap
        
          
    def final_energy(self):
        """Final total energy of the calculation read from vasprun.xml with Pymatgen"""
        final_energy = None
        if self.is_converged:
            if 'Vasprun' in self.outputs.keys():
                if self.vasprun:
                    final_energy = self.vasprun.final_energy
        return final_energy
     
            
    def final_structure(self):
        """Final structure read from "vasprun.xml" with Pymatgen"""
        if self.vasprun:
            final_structure = self.vasprun.structures[-1]
        else:
            final_structure = None
        return final_structure
                    

    def get_inputs(self,sync=False):
        """
        Read VaspInput from directory
        """
        if sync:
            self.sync_from_hpc()
        inputs = VaspInput.from_directory(self.path)
        self.inputs = inputs
        return
    
    
    def get_outputs(self,sync=False):
        """
        Get outputs dictionary from the data stored in the job directory. "vasprun.xml" is 
        read with Pymatgen
        """
        if sync:
            self.sync_from_hpc()
        path = self.path
        outputs = {}
        if op.isfile(op.join(path,'vasprun.xml')):
            try:
                outputs['Vasprun'] = Vasprun(op.join(path,'vasprun.xml'))
            except:
                print('Warning: Reading of vasprun.xml in "%s" failed'%path)
                outputs['Vasprun'] = None
        self.outputs = outputs
        return


    def hubbard(self):
        """
        Generate dictionary with U paramenters from LDAUU tag in INCAR file

        Returns
        -------
        U_dict : (dict)
            Dictionary with Elements as keys and U parameters as values.
        """
        U_dict = {}
        incar = self.incar
        if 'LDAUU' in incar.keys():
            ldauu = incar['LDAUU']
            elements = self.initial_structure.composition.elements
            if isinstance(ldauu,str):
                ldauu = ldauu.split()
            for i in range(0,len(ldauu)):
                U_dict[elements[i]] = ldauu[i]
        else:
            print('No LDAUU tag present in INCAR in Job "%s"' %self.name)
            
        return U_dict
            

    def nelectrons(self):
        """
        Number of electrons in the system. If 'NELECT' tag is in INCAR that value is returned.
        Else the sum of valence electrons from POTCAR is returned.
        """
        if 'NELECT' in self.incar.keys():
            nelect = self.incar['NELECT']
        else:
            val = {}
            for p in self.potcar:
                val[p.element] = p.nelectrons
            nelect = sum([ val[el.symbol]*self.initial_structure.composition[el] 
                           for el in self.initial_structure.composition])
        return nelect

    
    def write_input(self):
        """Write "inputs" dictionary to files. The VaspInput class from Pymatgen is used."""
        script_handler = ScriptHandler(**self.job_settings)
        script_handler.write_script(path=self.path)
        inputs = self.inputs
        inputs.write_input(output_dir=self.path,make_dir_if_not_present=True)
        return
    
    

    
class VaspNEBJob(Job):
    
    
    @staticmethod
    def from_directory(path,job_script_filename='job.sh'):
        """
        Builds VaspNEBjob object from data stored in a directory. Inputs dict is constructed
        by reading with Pymatgen INCAR, KPOINTS and POTCAR and creating a series of Structure 
        objects read from POSCARs in the images folders. 
        Inputs is thus a dict with "structures", "INCAR","KPOINTS","POTCAR" as keys.
        Output files are read usign Pymatgen NEBAnalysis and Vasprun classes.
        Job settings are read from the job script file.

        Parameters
        ----------
        path : (str)
            Path were job data is stored.
        job_script_filename : (str), optional
            Filename of job script. The default is 'job.sh'.

        Returns
        -------
        VaspNEBJob object.
        
        """                
        inputs = {}
        structures = []
        path = op.abspath(path)
        dirs = [d[0] for d in os.walk(path)]
        for d in dirs:
            image_name = op.relpath(d,start=path)
            if all(c.isdigit() for c in list(image_name)): #check if folder is image (all characters in folder rel path need to be numbers)
                image_path = d
                structure = Poscar.from_file(op.join(image_path,'POSCAR')).structure
                structures.append(structure)

        inputs['structures'] = structures           
        inputs['INCAR'] = Incar.from_file(op.join(path,'INCAR'))
        inputs['KPOINTS'] = Kpoints.from_file(op.join(path,'KPOINTS'))
        inputs['POTCAR'] = Potcar.from_file(op.join(path,'POTCAR'))
        
        outputs = {}
        if op.isfile(op.join(path,'vasprun.xml')):
            try:
                outputs['Vasprun'] = Vasprun(op.join(path,'vasprun.xml'))
            except:
                print('Warning: Reading of vasprun.xml in "%s" failed'%path)
                outputs['Vasprun'] = None   
        try:
            outputs['NEBAnalysis'] = NEBAnalysis.from_dir(path)
        except:
            print('Warning: NEB output reading with NEBAnalysis in "%s" failed'%path)
            outputs['NEBAnalysis'] = None
            
        s = ScriptHandler.from_file(path,filename=job_script_filename)
        job_settings =  s.settings
        
        return VaspNEBJob(path,inputs,job_settings,outputs)

 
    def delete_outputs(self,safety=True):
        """
        Delete files that aren't input files (INCAR,KPOINTS,POSCAR,POTCAR)
        """
        if safety:
            inp = input('Are you sure you want to delete outputs of Job %s ?: (y/n)' %self.name)
            if inp in ('y','Y'):
                delete = True
            else:
                delete = False
        else:
            delete= True
            
        if delete:
            dirs = self.image_dirs
            dirs.append(self.path)
            for d in dirs:                
                files = [f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]
                for f in files:
                    if f not in ['INCAR','KPOINTS','POSCAR','POTCAR',self.job_script_filename]:
                        os.remove(os.path.join(d,f))
                        print('Deleted file %s'%os.path.join(d,f))   
        return
    
    @property
    def images(self):
        return len(self.inputs['structures'])-2
    
    @property
    def image_dirs(self):
        """
        Directories of images for NEB calculations. Directories are selected if all characters in the
        directory name are digits.
        """
        dirs = []
        path = self.path
        path = op.abspath(path)
        for d in os.walk(path):
            directory = d[0]
            image_name = op.relpath(directory,start=path)
            if all(c.isdigit() for c in list(image_name)): #check if folder is image (all characters in folder rel path need to be numbers)
                dirs.append(directory)
        dirs.sort()
        return dirs
    
    @property
    def structures(self):
        return self.inputs['structures']

    @property
    def incar(self):
        return self.inputs['INCAR']
    
    @property
    def kpoints(self):
        return self.inputs['KPOINTS']
    
    @property
    def potcar(self):
        return self.inputs['POTCAR']
    
    @property
    def vasprun(self):
        if 'Vasprun' in self.outputs.keys():
            return self.outputs['Vasprun']
        else:
            raise ValueError('"vasprun.xml" file is not present in Job directory')
    
    @property
    def neb_analysis(self):
        return self.outputs['NEBAnalysis']
    

    @property
    def formula(self):
        """Complete formula from initial structure (read with Pymatgen)"""
        if self.initial_structure:
            return self.initial_structure.composition.formula
        else:
            return None

        
    @property
    def initial_structure(self):
        """Initial structure read from first element of ""structures" attribute. """
        return self.structures[0]

    
    @property
    def is_converged(self):
        """
        Reads Pymatgen Vasprun object and returns "True" if the calculation is converged,
        "False" if reading failed, and "None" if is not present in the outputs dictionary.
        """
        is_converged = None
        if self.outputs:
            if 'Vasprun' in self.outputs.keys():
                is_converged = False
                if self.vasprun:
                    vasprun = self.vasprun
                    conv_el, conv_ionic = False, False
                    if vasprun:
                        conv_el = vasprun.converged_electronic
                        conv_ionic = vasprun.converged_ionic
                    if conv_el and conv_ionic:
                        is_converged = True                
        if not is_converged:
            if self.is_step_limit_reached:
                is_converged = True
            
        return is_converged    
    
    
    @property
    def is_step_limit_reached(self):
        """
        Reads number of ionic steps from the OSZICAR file with Pymatgen and returns True if 
        is equal to the step limit in INCAR file (NSW tag)
        """
        limit_reached = True
        image_dirs = self.image_dirs
        for d in image_dirs:
            if d != image_dirs[0] and d != image_dirs[-1]:
                if not os.path.isfile(os.path.join(d,'OSZICAR')): # check if OSZICAR files are present 
                    limit_reached = False
                else:                
                    n_steps = len(Oszicar(os.path.join(d,'OSZICAR')).ionic_steps)
                    nsw = Incar.from_file(op.join(self.path,'INCAR'))['NSW'] # check NSW from INCAR in parent directory
                    if nsw != n_steps:
                        limit_reached = False
        return limit_reached
    

    def charge(self):
        """
        Charge of the system calculated as the difference between the value of "NELECT"
        in the INCAR and the number of electrons in POTCAR. If "NELECT" is not present 
        charge is set to 0.
        """
        charge = 0
        if 'NELECT' in self.incar.keys():
            nelect = self.incar['NELECT']
            val = {}
            for p in self.potcar:
                val[p.element] = p.nelectrons
            neutral = sum([ val[el.symbol]*self.initial_structure.composition[el] 
                           for el in self.initial_structure.composition])
            charge = neutral - nelect
        return charge


    def get_inputs(self,sync=False):
        """
        Read inputs from Job directory
        """
        if sync:
            self.sync_from_hpc()
        inputs = {}
        structures = []
        path = op.abspath(self.path)
        dirs = [d[0] for d in os.walk(path)]
        for d in dirs:
            image_name = op.relpath(d,start=path)
            if all(c.isdigit() for c in list(image_name)): #check if folder is image (all characters in folder rel path need to be numbers)
                image_path = d
                structure = Poscar.from_file(op.join(image_path,'POSCAR')).structure
                structures.append(structure)

        inputs['structures'] = structures           
        inputs['INCAR'] = Incar.from_file(op.join(path,'INCAR'))
        inputs['KPOINTS'] = Kpoints.from_file(op.join(path,'KPOINTS'))
        inputs['POTCAR'] = Potcar.from_file(op.join(path,'POTCAR'))
        
        self.inputs = inputs
        return


    def get_outputs(self,sync=False):
        """
        Read outputs from Job directory
        """
        if sync:
            self.sync_from_hpc()
        outputs = {}
        path = self.path
        if op.isfile(op.join(path,'vasprun.xml')):
            try:
                outputs['Vasprun'] = Vasprun(op.join(path,'vasprun.xml'))
            except:
                print('Warning: Reading of vasprun.xml in "%s" failed'%path)
                outputs['Vasprun'] = None   
        try:
            outputs['NEBAnalysis'] = NEBAnalysis.from_dir(path)
        except:
            print('Warning: NEB output reading with NEBAnalysis in "%s" failed'%path)
            outputs['NEBAnalysis'] = None
            
        self.outputs = outputs
        return
    

    def nelectrons(self):
        """
        Number of electrons in the system. If 'NELECT' tag is in INCAR that value is returned.
        Else the sum of valence electrons from POTCAR is returned.
        """
        if 'NELECT' in self.incar.keys():
            nelect = self.incar['NELECT']
        else:
            val = {}
            for p in self.potcar:
                val[p.element] = p.nelectrons
            nelect = sum([ val[el.symbol]*self.initial_structure.composition[el] 
                           for el in self.initial_structure.composition])
        return nelect

    
    def write_input(self,write_structures=True):
        """
        Write input files in all image directories
        """
        path = op.abspath(self.path)
        
        self.job_settings['nodes'] = self.images
        if 'add_automation' not in self.job_settings:
            self.job_settings['add_automation'] = None
               
        incar = self.inputs['INCAR']
        kpoints = self.inputs['KPOINTS']
        potcar = self.inputs['POTCAR']
        job_settings = self.job_settings

        if write_structures:
            self.write_structures()
        
        incar.write_file(op.join(path,'INCAR'))
        kpoints.write_file(op.join(path,'KPOINTS'))
        potcar.write_file(op.join(path,'POTCAR'))
        ScriptHandler(**job_settings).write_script(path=path)

    
    def write_structures(self):
        """
        Writes POSCAR files in image directories
        """
        path = self.path
        structures = self.inputs['structures']
        for s in structures:
            index = structures.index(s)
            image_path = op.join(path,str(index).zfill(2)) #folders will be named 00,01,..,XX
            if not op.exists(image_path):
                os.makedirs(image_path)
            Poscar(s).write_file(op.join(image_path,'POSCAR'))
        return