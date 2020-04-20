
import os
import os.path as op
from my_functions.data.jobs import VaspJob
import pandas as pd


class Dataset:
    
    def __init__(self,path=None,name=None,jobs=None): 
        
        self.path = path if path else os.getcwd()
        self.name = name if name else op.basename(op.abspath(self.path))
        self.jobs = jobs
        self.groupped_jobs = self.group_jobs()

    @staticmethod
    def from_directory(path,job_script_filename='job.sh'): #to change job script filename in 'job.sh'
        jobs = []
        for root , dirs, files in os.walk(path):
            if files != [] and job_script_filename in files:
                if ('INCAR' and 'KPOINTS' and 'POSCAR' and 'POTCAR') in files:
                    j = VaspJob.from_directory(root)
                    j.job_script_filename = job_script_filename
                    jobs.append(j)
     
        return  Dataset(path=path,jobs=jobs)
        
    @property
    def groups(self):
        return next(os.walk(self.path))[1]
    
    def group_jobs(self):
        path = op.abspath(self.path)
        gjobs = {}
        groups = self.groups
        for group in groups:
            gjobs[group] = {}
            group_path = op.join(path,group)
            for job in self.jobs:
                if group in job.path:
                    jpath = op.abspath(job.path)
                    nodes = jpath.replace(op.commonpath([group_path,jpath]),'')
                    gjobs[group][nodes] = job
                    job.group = group
                    job.nodes = nodes                    
        return gjobs
    
    
    def jobs_table(self,jobs=[],properties_to_display=[],common_node=None):
        
        jobs = jobs if jobs else self.jobs
        table = []
        index = []
        for j in jobs:
            d = {}
            d['formula'] = j.formula
            d['group'] = j.group
            d['nodes'] = j.nodes
            d['is_converged'] = j.is_converged
            for feature in properties_to_display:
                d[feature] = getattr(j,feature) ()
            if common_node:
                if common_node in j.nodes:
                    table.append(d)
                    index.append(j.name)
            else:
                table.append(d)
                index.append(j.name)
            
        df = pd.DataFrame(table,index=index)
        df.index.name = 'job_name'
        
        return df
            
            
            

            
                
        
                

            
        
        