#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:33:19 2023

@author: villa
"""
from pynter.testing.core import PynterTest
from pynter.testing.structure import StructureTest
from pynter.testing.hpc import JobSettingsTest
from pynter.testing.vasp import VaspInputsTest, VaspOutputsTest

class JobTest(PynterTest):
    """
    Contains methods to test Job objects
    """
    
    def __init__(self,jobclass='VaspJob'):
        """
        Specify which class the job object belongs to. Currently only "VaspJob" is available.
        """
        if jobclass == 'VaspJob':
            self.cls = VaspJobTest
        elif jobclass == 'VaspNEBJob':
            self.cls = VaspNEBJobTest
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
        #self.assert_settings_equal(job1,job2)  # don't test for equal settings since they are system dependent
            
            
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
            
            
class VaspJobTest(JobTest):
    """
    Provides methods to test VaspJob objects
    """
    def assert_inputs_equal(self,job1,job2,include_incar=True,system_only=False):
        VaspInputsTest().assert_VaspInput_equal(job1.inputs,job2.inputs,include_incar,system_only)

    def assert_outputs_equal(self,job1,job2):
        VaspOutputsTest().assert_ComputedEntry_equal(job1.computed_entry,job2.computed_entry)
        if "Vasprun" in [job1.outputs.keys(),job2.outputs.keys()]:
            VaspOutputsTest().assert_Vasprun_equal(job1.vasprun, job2.vasprun)
            
            
class VaspNEBJobTest(JobTest):
    """
    Provides methods to test VaspNEBJob objects
    """
    def assert_inputs_equal(self,job1,job2,include_incar=True,system_only=False):
        """
        Assert incar, kpoints, potcar and images.
        """
        VaspInputsTest().assert_Incar_equal(job1.incar, job2.incar, system_only=system_only)
        VaspInputsTest().assert_Kpoints_equal(job1.kpoints, job2.kpoints)
        VaspInputsTest().assert_Potcar_equal(job1.potcar, job2.potcar)

        if len(job1.structures) != len(job2.structures):
            raise AssertionError('Number of images differ')
        else:
            for i in range(0,len(job1.structures)):
                StructureTest().assert_Structure_equal(
                    job1.structures[i], job2.structures[i])
            
    def assert_outputs_equal(self,job1,job2):
        """
        Assert pymatgen NEBAnalysis object
        """
        self.assert_object_almost_equal(job1.neb_analysis.as_dict(), job2.neb_analysis.as_dict())           
            
            
            

