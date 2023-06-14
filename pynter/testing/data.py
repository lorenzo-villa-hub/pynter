#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""


import importlib

from pynter.testing.core import PynterTest
from pynter.testing.slurm import JobSettingsTest

class JobTest(PynterTest):
    """
    Contains methods to test Job objects
    """
    
    def __init__(self,jobclass='VaspJob'):
        """
        Specify which class the job object belongs to. Currently only "VaspJob" is available.
        """
        if jobclass == 'VaspJob':
            module = importlib.import_module("pynter.testing.vasp")
            self.cls = getattr(module,'VaspJobTest')
        elif jobclass == 'VaspNEBJob':
            module = importlib.import_module("pynter.testing.vasp")
            self.cls = getattr(module,'VaspNEBJobTest')
        super().__init__()
        
    def assert_inputs_equal(self,job1,job2,**kwargs):
        self.cls().assert_inputs_equal(job1,job2,**kwargs)
        return

    def assert_outputs_equal(self,job1,job2,**kwargs):
        self.cls().assert_outputs_equal(job1,job2,**kwargs)
        return
        
    def assert_settings_equal(self,job1,job2):
        JobSettingsTest().assert_job_settings_equal(job1.job_settings,job2.job_settings)
        
    def assert_job_equal(self,job1,job2,kwargs_input={},kwargs_output={}):
        self.assert_inputs_equal(job1,job2,**kwargs_input)
        if job1.outputs or job2.outputs:
            self.assert_outputs_equal(job1,job2,**kwargs_output)
        self.assert_settings_equal(job1,job2)
            
            
class DatasetTest(PynterTest):
    
    def assert_dataset_equal(self,ds1,ds2,kwargs_input={},kwargs_output={}):
        """
        Check if jobs in dataset are equal. Jobs are sorted before they are checked.
        """
        ds1.sort_jobs()
        ds2.sort_jobs()
        for i in range(len(ds1.jobs)):
            job1,job2 = ds1.jobs[i], ds2.jobs[i]
            JobTest(jobclass=job1.jobclass).assert_job_equal(job1, job2,kwargs_input,kwargs_output)
            
    def assert_jobs_equal(self,jobs1,jobs2,kwargs_input={},kwargs_output={}):
        """
        Check if lists of jobs are equal.
        """
        assert len(jobs1) == len(jobs2), "Job lists have different lengths"
        for i in range(0,len(jobs1)):
            j1, j2 = jobs1[i], jobs2[i]
            assert j1.jobclass == j2.jobclass, "Jobs belong to different classes" 
            JobTest(jobclass=j1.jobclass).assert_job_equal(j1, j2,kwargs_input,kwargs_output)
            
            
            
            
            
            

