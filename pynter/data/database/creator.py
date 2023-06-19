#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 11:12:09 2020

@author: villa
"""

import warnings
from matgendb.creator import VaspToDbTaskDrone
from pynter import SETTINGS

dbconfig = SETTINGS['dbconfig']

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
    def __init__(self, job, **kwargs):
        
        for k,v in dbconfig['vasp'].items():
            if k not in kwargs:
                kwargs[k] = v
                   
        super().__init__(use_full_uri=False,**kwargs)

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
                warnings.warn(f'VaspJob "{self.job.name}" is not converged and will not be added into the database',UserWarning)
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
    
        
        