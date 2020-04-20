
import os
import os.path as op
from pymatgen.io.vasp.inputs import VaspInput, Poscar
from pymatgen.io.vasp.outputs import Vasprun


class Job:
    
    def __init__(self,path=None,inputs=None,job_settings=None,outputs=None):
        
        self.path = path if path else os.getcwd()
        self.inputs = inputs
        self.job_settings = job_settings
        self.outputs = outputs

    @property 
    def name(self): # here use from_file method of ScriptHandler in job_script module
         return 'test' # temporary name to test 

     
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
    def formula(self):
        return self.initial_structure.composition.formula
        
    @property
    def initial_structure(self):
        return Poscar.from_file(op.join(self.path,'POSCAR')).structure
    
    
    def charge(self):
        nelect = self.inputs['INCAR']['NELECT']
        val = {}
        for p in self.inputs['POTCAR']:
            val[p.element] = p.nelectrons
        neutral = sum([ val[el.symbol]*self.initial_structure.composition[el] 
                       for el in self.initial_structure.composition])
        return neutral - nelect
    
      
    def final_energy(self):
        if self.outputs['vasprun']:
            final_energy = self.outputs['vasprun'].final_energy
        else:
            final_energy = None
        return final_energy
     
            
    def final_structure(self):
        if self.outputs['vasprun']:
            final_structure = self.outputs['vasprun'].structures[-1]
        else:
            final_structure = None
        return final_structure
    

        
    