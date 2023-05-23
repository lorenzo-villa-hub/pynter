#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""

import unittest
import importlib


class CompareJobs(unittest.TestCase):
    
    def __init__(self,jobclass='VaspJob'):
        if jobclass == 'VaspJob':
            module = importlib.import_module("pynter.vasp.tests.compare")
            self.cls = getattr(module,'CompareVaspJobs')
        super().__init__()
        
    def compare_inputs(self,job1,job2,**kwargs):
        self.cls().compare_inputs(job1,job2,**kwargs)
        return

    def compare_outputs(self,job1,job2,**kwargs):
        self.cls().compare_outputs(job1,job2,**kwargs)
        return
        
    def compare_settings(self,job1,job2):
        self.assertEqual(job1.job_settings,job2.job_settings)
        self.assertEqual(job1.name,job2.name)
        
    def compare(self,job1,job2,kwargs_input={},kwargs_output={}):
        self.compare_inputs(job1,job2,**kwargs_input)
        if job1.outputs or job2.outputs:
            self.compare_outputs(job1,job2)
        self.compare_settings(job1,job2)
            
            
class CompareDatasets(unittest.TestCase):
    
    def compare_jobs(self,ds1,ds2,kwargs_input={},kwargs_output={}):
        ds1.sort_jobs()
        ds2.sort_jobs()
        for i in range(len(ds1.jobs)):
            job1,job2 = ds1.jobs[i], ds2.jobs[i]
            CompareJobs(jobclass=job1.jobclass).compare(job1, job2,kwargs_input,kwargs_output)
            
            

