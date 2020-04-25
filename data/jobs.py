
from abc import ABCMeta, abstractmethod
import os
import os.path as op
from pymatgen.io.vasp.inputs import VaspInput, Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pynter.slurm.job_script import ScriptHandler
from pynter.slurm.interface import HPCInterface
from pynter.tools.grep import grep_list
from pynter.__init__ import load_config


class Job:
    
    def __init__(self,path=None,inputs=None,job_settings=None,outputs=None,job_script_filename='job.sh'):
        
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
        
        
    @property 
    def name(self): 
        s = ScriptHandler.from_file(self.path,filename=self.job_script_filename)
        return s.settings['name'] 

        
    def cancel_job(self):
        
        hpc = HPCInterface()
        job_id = self.job_id()
        hpc.cancel_jobs(job_id)
        
        return 


    @abstractmethod
    def get_outputs(self):
        pass


    def job_id(self):
        
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
        
        hpc = HPCInterface()
        stdout,stderr = hpc.qstat()
        
        return stdout,stderr
        

    def run_job(self,write_input=True):
        
        if write_input:
            self.write_input()
        hpc = HPCInterface()
        hpc.rsync_to_hpc()
        stdout,stderr = hpc.sbatch(path=self.path_in_hpc,job_script_filename=self.job_script_filename)
        
        return stdout,stderr
    

    def sync_job(self,stdouts=False):
        hpc = HPCInterface()
        abs_path = op.abspath(self.path)
        localdir = abs_path #abs_path.replace(op.basename(abs_path),'')
        stdout,stderr = hpc.rsync_from_hpc(remotedir=self.path_in_hpc,localdir=localdir)
        if stdouts:
            return stdout,stderr
        else:
            return
        


    def status(self):
        
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
            
        return status
 
           
    @abstractmethod
    def write_input():
        pass
        
     
class VaspJob(Job):
 
      
    @staticmethod
    def from_directory(path,job_script_filename='job.sh'):
                
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
        return self.initial_structure.composition.formula

        
    @property
    def initial_structure(self):
        return Poscar.from_file(op.join(self.path,'POSCAR')).structure
    
    
    @property
    def is_converged(self):
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
 
    
    def get_outputs(self):
        
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
        final_energy = None
        if 'vasprun' in self.outputs.keys():
            if self.outputs['vasprun']:
                final_energy = self.outputs['vasprun'].final_energy
        return final_energy
     
            
    def final_structure(self):
        if self.outputs['vasprun']:
            final_structure = self.outputs['vasprun'].structures[-1]
        else:
            final_structure = None
        return final_structure
                    
    
    def write_input(self):
        
        script_handler = ScriptHandler(**self.job_settings)
        script_handler.write_script(path=self.path)
        inputs = self.inputs
        inputs.write_input(output_dir=self.path,make_dir_if_not_present=True)
        return