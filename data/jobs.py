
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
        
        self.path = path if path else os.getcwd()
        self.inputs = inputs
        self.job_settings = job_settings
        self.outputs = outputs
        self.job_script_filename = job_script_filename
        
        self.localdir = HPCInterface().localdir
        self.workdir = HPCInterface().workdir
        self.path_relative = op.abspath(self.path).replace(self.localdir,'')
        
        
    @property 
    def name(self): 
        s = ScriptHandler.from_file(self.path,filename=self.job_script_filename)
        return s.settings['name'] 

        

    def cancel_job(self):
        
        hpc = HPCInterface()
        job_id = self.job_id()
        hpc.cancel_jobs(job_id)
        
        return 


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
        

    def run_job(self):
        
        hpc = HPCInterface()
        hpc.rsync_to_hpc()
        hpc_path = self.workdir + self.path_relative
        stdout,stderr = hpc.sbatch(path=hpc_path,job_script_filename=self.job_script_filename)
        
        return stdout,stderr
    

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
                    