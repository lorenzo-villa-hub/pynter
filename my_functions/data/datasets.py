
import os
import os.path as op
from pymatgen.io.vasp.inputs import VaspInput, Poscar
from pymatgen.io.vasp.outputs import Vasprun

class Dataset:
    
    def __init__(self,path=None,name=None,jobs=None): 
        
        self.path = path if path else os.getcwd()
        self.name = name if name else op.basename(self.path)
        self.jobs = jobs

    @staticmethod
    def from_directory(path,job_script_filename='job_vasp.sh'): #to change job script filename in 'job.sh'
        jobs = []
        for root , dirs, files in os.walk(path):
            if files != [] and job_script_filename in files:
                if ('INCAR' and 'KPOINTS' and 'POSCAR' and 'POTCAR') in files:
                    jobs.append(VaspJob.from_directory(root))
     
        return  Dataset(path=path,jobs=jobs)
        
    @property
    def groups(self):
        return next(os.walk(self.path))[1]
    

                
class Job:
    
    def __init__(self,path=None,inputs=None,job_settings=None,outputs=None):
        
        self.path = path if path else os.getcwd()
        self.inputs = inputs
        self.job_settings = job_settings
        self.outputs = outputs


class VaspJob(Job):
       
    @staticmethod
    def from_directory(path):
                
        inputs = VaspInput.from_directory(path)
        outputs = {}
        if op.isfile(op.join(path,'vasprun.xml')):
            try:
                outputs['vasprun'] = Vasprun(op.join(path,'vasprun.xml'))
            except:
                print('Warning: Reading of vasprun.xml in "%s" failed'%path)
                outputs['vasprun'] = None
        
        job_settings = None # here use from_file for job settings (still to write)
        
        return VaspJob(path,inputs,job_settings,outputs)
            
    @property
    def final_energy(self):
        if self.outputs['vasprun']:
            final_energy = self.outputs['vasprun'].final_energy
        else:
            final_energy = None
        return final_energy
     
    @property 
    def initial_structure(self):
        return Poscar.from_file(op.join(self.path,'POSCAR')).structure
    
    @property
    def final_structure(self):
        if self.outputs['vasprun']:
            final_structure = self.outputs['vasprun'].structures[-1]
        else:
            final_structure = None
        return final_structure
            
        
        