#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:55:50 2020

@author: villa
"""


from matgendb.query_engine import QueryEngine
from pymatgen.io.vasp.inputs import VaspInput, Incar,Kpoints,Poscar,Potcar
from pymatgen.core.structure import Structure
from pynter.data.jobs import VaspJob


class VaspJobQuery(QueryEngine):
    
    def __init__(self, host="127.0.0.1", port=27017, database="vasp",
                 user=None, password=None, collection="tasks",
                 aliases_config=None, default_properties=None,
                 query_post=None, result_post=None,
                 connection=None, replicaset=None, **ignore):
        
        super().__init__(host,port,database,user,password,collection,aliases_config,
                         default_properties,query_post,result_post,connection,replicaset,**ignore)
        
        self.optional_data = ['calculations','dir_name','final_energy','eigenvalue_band_properties',
                         'job_settings','job_script_filename','job_name','is_converged','structures']
        

    def convert_computed_entry_to_job(self,entry):
        e = entry
        path = e.data['dir_name']
        
        inp = e.data['calculations'][0]['input']           
        incar = Incar(inp['incar'])
        kpoints = Kpoints.from_dict(inp['kpoints'])
        poscar = Poscar(e.structure)
        potcar = Potcar(inp['potcar'])
        
        inputs = VaspInput(incar, kpoints, poscar, potcar)
        job_settings = e.data['job_settings']
        job_script_filename = e.data['job_script_filename']
        name = e.data['job_name']
        
        outputs = {'ComputedStructureEntry':e}
        
        vaspjob =  VaspJob(path,inputs,job_settings,outputs,job_script_filename,name)
        vaspjob._is_converged = e.data['is_converged']
        vaspjob._band_structure = None
            
        return vaspjob            


    def get_entries_optional_data(self,criteria,optional_data=[]):

        optional_data = optional_data + self.optional_data
        entries = self.get_entries(criteria,inc_structure=True,optional_data=optional_data)
        return entries


    def get_entries_in_system_optional_data(self,elements,optional_data=[],additional_criteria=None):
        
        optional_data = optional_data + self.optional_data
        entries = self.get_entries_in_system(elements,inc_structure=True, optional_data=optional_data,
                                             additional_criteria=additional_criteria)
        return entries
    
    
    def get_jobs(self,criteria,optional_data=[]):
        
        optional_data = optional_data + self.optional_data
        jobs = []
        entries = self.get_entries_optional_data(criteria,optional_data=optional_data)
        for e in entries:
            j = self.convert_computed_entry_to_job(e)
            jobs.append(j)
        
        return jobs


            