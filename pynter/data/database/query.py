#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:55:50 2020

@author: villa
"""


from matgendb.query_engine import QueryEngine
from pymatgen.io.vasp.inputs import VaspInput, Incar,Kpoints,Poscar,Potcar
from pynter.vasp.jobs import VaspJob
from pynter import SETTINGS


dbconfig = SETTINGS['dbconfig']

class VaspJobQuery(QueryEngine):
    """
    Subclass of QueryEngine in pymatgen-db (matgendb) package.
    To set up the database for VaspJob we use the pymatgen-db package which already provides a useful
    interface to the other pymatgen objects that are already used (like Structure or ComputedEntry).
    This class helps to interface a VaspJob object with the database.
    This class sets up the default "optional_data" to get from the db which allows to convert 
    pymatgen's ComputedStructureEntry objects retrieved from the db to VaspJob objects. 
    """
    def __init__(self,aliases_config=None,default_properties=None,query_post=None,
                 result_post=None,connection=None,replicaset=None,**kwargs):
        
        for k,v in dbconfig['vasp'].items():
            if k not in kwargs:
                kwargs[k] = v
        host = kwargs['host']
        port = kwargs['port']
        database = kwargs['database']
        user = kwargs['user']
        password = kwargs['password']
        collection = kwargs['collection']
        
        #can't pass kwargs directly because of **ignore in QueryEngine
        super().__init__(host=host,port=port,database=database,user=user,password=password,collection=collection,
                         aliases_config=aliases_config,default_properties=default_properties,query_post=query_post,
                         result_post=result_post,connection=connection,replicaset=replicaset) 
        
        self.optional_data = ['calculations','dir_name','final_energy','eigenvalue_band_properties',
                         'job_settings','job_script_filename','job_name','is_converged','structures']
        

    def convert_computed_entry_to_job(self,entry):
        """
        Convert ComputedStructureEntry into VaspJob object

        Parameters
        ----------
        entry : 
            ComputedStructureEntry.

        Returns
        -------
        vaspjob : 
            VaspJob object.
        """
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
        """
        Get entries with get_entries method with default "optional_data".

        Parameters
        ----------
        criteria : 
            Criteria obeying the same syntax as query.
        optional_data : (list), optional
            Optional data to include in the entry. The default is [].

        Returns
        -------
        entries : (list)
            List of ComputedStructureEntry objects.
        """
        optional_data = optional_data + self.optional_data
        entries = self.get_entries(criteria,inc_structure=True,optional_data=optional_data)
        return entries


    def get_entries_in_system_optional_data(self,elements,optional_data=[],additional_criteria=None):
        """
        Uses get_entries_in_system method with default "optional_data".
        Gets all entries in a chemical system, e.g. Li-Fe-O will return all
        Li-O, Fe-O, Li-Fe, Li-Fe-O compounds.
        
        Parameters
        ----------
        elements: (list)
            List of strings with element symbols.
        optional_data : (list), optional
            Optional data to include in the entry. The default is [].
        criteria : 
            Added ability to provide additional criteria other than just
            the chemical system.

        Returns
        -------
        entries : (list)
            List of ComputedStructureEntry objects.
        """        
        optional_data = optional_data + self.optional_data
        entries = self.get_entries_in_system(elements,inc_structure=True, optional_data=optional_data,
                                             additional_criteria=additional_criteria)
        return entries
    
    
    def get_jobs(self,criteria,optional_data=[]):
        """
        Get VaspJob objects from db. First retrives ComputedStructureEntry objects than 
        converts them into VaspJob objects.

        Parameters
        ----------
        criteria : 
            Criteria obeying the same syntax as query.
        optional_data : (list), optional
            Optional data to include in the entry. The default is [].

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        optional_data = optional_data + self.optional_data
        jobs = []
        entries = self.get_entries_optional_data(criteria,optional_data=optional_data)
        for e in entries:
            j = self.convert_computed_entry_to_job(e)
            jobs.append(j)
        
        return jobs


            