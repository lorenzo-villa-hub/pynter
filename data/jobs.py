
from abc import ABCMeta, abstractmethod
import os
import os.path as op
from pymatgen.io.vasp.inputs import VaspInput, Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter,BSDOSPlotter
from pynter.slurm.job_script import ScriptHandler
from pynter.slurm.interface import HPCInterface
from pynter.tools.grep import grep_list
from pynter.__init__ import load_config


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
        
        __metaclass__ = ABCMeta
        
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
        else:
            s = ScriptHandler.from_file(self.path,filename=self.job_script_filename)
            self.name = s.settings['name']
        if not self.job_settings:
            self.job_settings = {}
        self.job_settings['name'] = self.name


    def __str__(self):
        printout = 'Job "%s" of group "%s"' %(self.name, self.group)
        return printout
    
    def __repr__(self):
        return self.__str__()
        
        
    def cancel_job(self):
        """Cancel job on HPC"""
        hpc = HPCInterface()
        job_id = self.job_id()
        hpc.cancel_jobs(job_id)
        
        return 


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
                outputs['vasprun'] = Vasprun(op.join(path,'vasprun.xml'))
            except:
                print('Warning: Reading of vasprun.xml in "%s" failed'%path)
                outputs['vasprun'] = None
        
        s = ScriptHandler.from_file(path,filename=job_script_filename)
        job_settings =  s.settings
        
        return VaspJob(path,inputs,job_settings,outputs)

    
    @property
    def formula(self):
        """Complete formula from initial structure (read with Pymatgen)"""
        if self.initial_structure:
            return self.initial_structure.composition.formula
        else:
            return None

        
    @property
    def initial_structure(self):
        """Initial structure read from Poscar file with Pymatgen"""
        poscar_file = op.join(self.path,'POSCAR')
        if op.isfile(poscar_file):
            return Poscar.from_file(poscar_file).structure
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
            if 'vasprun' in self.outputs.keys():
                is_converged = False
                if self.outputs['vasprun']:
                    vasprun = self.outputs['vasprun']
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
        if 'NELECT' in self.inputs['INCAR'].keys():
            nelect = self.inputs['INCAR']['NELECT']
            val = {}
            for p in self.inputs['POTCAR']:
                val[p.element] = p.nelectrons
            neutral = sum([ val[el.symbol]*self.initial_structure.composition[el] 
                           for el in self.initial_structure.composition])
            charge = neutral - nelect
        return charge
 
    
    def get_outputs(self,sync=False):
        """
        Get outputs dictionary from the data stored in the job directory. "vasprun.xml" is 
        read with Pymatgen
        """
        if sync:
            self.sync_job()
        path = self.path
        outputs = {}
        if op.isfile(op.join(path,'vasprun.xml')):
            try:
                outputs['vasprun'] = Vasprun(op.join(path,'vasprun.xml'))
            except:
                print('Warning: Reading of vasprun.xml in "%s" failed'%path)
                outputs['vasprun'] = None
        self.outputs = outputs
        return
        
          
    def final_energy(self):
        """Final total energy of the calculation read from vasprun.xml with Pymatgen"""
        final_energy = None
        if self.is_converged:
            if 'vasprun' in self.outputs.keys():
                if self.outputs['vasprun']:
                    final_energy = self.outputs['vasprun'].final_energy
        return final_energy
     
            
    def final_structure(self):
        """Final structure read from "vasprun.xml" with Pymatgen"""
        if self.outputs['vasprun']:
            final_structure = self.outputs['vasprun'].structures[-1]
        else:
            final_structure = None
        return final_structure
                    

    def plot_dos_bs(self):
        """
        Plot DOS and BS from data in vasprun.xml with Pymatgen
        """
        wdir = os.getcwd()
        os.chdir(self.path)
        if self.is_converged:
            vasprun = self.outputs['vasprun']
            bs = vasprun.get_band_structure(line_mode=True)
            dos = vasprun.complete_dos
            plt = BSDOSPlotter(bs_projection=None,dos_projection=None).get_plot(bs,dos)           
        else:
            raise ValueError(f'Job %s is not converged' %self.name)
        os.chdir(wdir)
        return plt
            
    
    def write_input(self):
        """Write "inputs" dictionary to files. the VaspInput class from Pymatgen is used."""
        script_handler = ScriptHandler(**self.job_settings)
        script_handler.write_script(path=self.path)
        inputs = self.inputs
        inputs.write_input(output_dir=self.path,make_dir_if_not_present=True)
        return