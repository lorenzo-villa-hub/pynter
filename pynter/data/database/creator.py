#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 11:12:09 2020

@author: villa
"""


from matgendb.creator import VaspToDbTaskDrone
from pynter.data.jobs import VaspJob
from pynter.data.datasets import Dataset

# for VaspJob we use pymatgen-db interface


class VaspJobDrone(VaspToDbTaskDrone):
    """
    Subclass of VaspToDbTaskDrone in pymatgen-db (matgendb) package.
    To set up the database for VaspJob we use the pymatgen-db package which already provides a useful
    interface to the other pymatgen objects that are already used (like Structure or ComputedEntry).
    This class helps to interface a VaspJob object with the database.
    For documentation on input args (https://github.com/materialsproject/pymatgen-db/blob/master/matgendb/creator.py).
    Warning: 'use_full_uri' has been set to False here to allow interface with VaspJob objects.
    """
    
    def __init__(self, job, host="127.0.0.1", port=27017, database="vasp",
                 user=None, password=None, collection="tasks",
                 parse_dos=False, compress_dos=False,parse_projected_eigen=False,
                 simulate_mode=False, additional_fields=None, update_duplicates=True,
                 mapi_key=None, use_full_uri=False, runs=None):
        
        super().__init__(host,port,database,user,password,collection,parse_dos,compress_dos,
                         parse_projected_eigen,simulate_mode,additional_fields,update_duplicates,
                         mapi_key,use_full_uri,runs)

        self.job = job
        self._path = self.job.path
    
    
    def assimilate_job(self,check_convergence=True):
        """
        Insert VaspJob intro database.
        """
        if check_convergence:
            if self.job.is_converged:
                doc = self.get_task_doc_from_files()     
                self._insert_doc(doc)
            else:
                print(f'VaspJob "{self.job.name}" is not converged and will not be added into the database')
        else:
            doc = self.get_task_doc_from_files()     
            self._insert_doc(doc)
        return
        
    
    def get_task_doc_from_files(self):
        """
        Get entire task doc from "vasprun.xml" in job path.
        """
        doc = self.get_task_doc(self._path)
        doc['job_settings'] = self.job.job_settings
        doc['job_script_filename'] = self.job.job_script_filename
        doc['job_name'] = self.job.name
        doc['is_converged'] = self.job.is_converged
        return doc
    
        
        