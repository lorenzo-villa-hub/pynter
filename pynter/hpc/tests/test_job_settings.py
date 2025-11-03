#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 13:46:42 2023

@author: villa
"""

import os

from pynter.hpc.slurm import Sbatch, JobSettings

from pynter.testing.core import PynterTest
from pynter.testing.hpc import JobSettingsTest


class TestJobSettings(PynterTest):
    
    def setUp(self):

        self.script_string = (
            "#!/bin/sh\n"
            "\n"
            "#SBATCH --account=projecttest0000\n"
            "#SBATCH --error=err.%j\n"
            "#SBATCH --job-name=test\n"
            "#SBATCH --mail-user=test@pynter.com\n"
            "#SBATCH --mem-per-cpu=3500\n"
            "#SBATCH --ntasks=96\n"
            "#SBATCH --output=out.%j\n"
            "#SBATCH --time=24:00:00\n"
            "\n"
            "module purge\n"
            "export APPLES\n"
            "export ORANGES\n"
            "ml intel/2020.4\n"
            "ml intelmpi/2020.4\n"
            "ml fftw/3.3.10\n"
            "test_HEADER\n"
            "test_HEADER2\n"
            "if [ ! -f POSCAR_initial ] ; then\n"
            "    cp POSCAR POSCAR_initial\n"
            "fi\n"
            "if [ -f CONTCAR ]\n"
            "then\n"
            "    cp CONTCAR POSCAR\n"
            "fi\n"
            "srun /home/test/code\n"
            "pynter analysis vasprun --convergence > convergence.txt\n"
            "if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then\n"
            "    automation.py\n"
            "    scancel ${SLURM_ARRAY_JOB_ID}_*\n"
            "fi\n"
            "test_BODY\n"
            "test_BODY2"
        )
        
        self.sbatch_kwargs = {
            'account': 'projecttest0000',
            'error': 'err.%j',
            'job-name': 'test',
            'mail-user': 'test@pynter.com',
            'mem-per-cpu':3500,
            'ntasks': 96,
            'output': 'out.%j', 
            'time': '24:00:00'
             }
        
        script_lines = [
            "module purge",
            "export APPLES",
            "export ORANGES",
            "ml intel/2020.4",
            "ml intelmpi/2020.4",
            "ml fftw/3.3.10",
            "test_HEADER",
            "test_HEADER2",
            "if [ ! -f POSCAR_initial ] ; then",
            "    cp POSCAR POSCAR_initial",
            "fi",
            "if [ -f CONTCAR ]",
            "then",
            "    cp CONTCAR POSCAR",
            "fi",
            "srun /home/test/code",
            "pynter analysis vasprun --convergence > convergence.txt",
            "if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then",
            "    automation.py",
            "    scancel ${SLURM_ARRAY_JOB_ID}_*",
            "fi",
            "test_BODY",
            "test_BODY2"            
            ]
        self.js = JobSettings(load_sbatch_defaults=False,sbatch=self.sbatch_kwargs,filename='job_test.sh',script_lines=script_lines)

        
    def test_string(self):
        assert self.js.get_bash_script() == self.script_string
        
    def test_from_file(self):
        js_from_file = JobSettings.from_bash_file(self.test_files_path,'job_test.sh')
        JobSettingsTest().assert_job_settings_equal(self.js,js_from_file)
        
    def test_write_script(self):
        self.js['filename'] = 'temp.sh'
        self.js.write_bash_file(self.test_files_path,filename='temp.sh')
        JobSettingsTest().assert_job_settings_equal(
            JobSettings.from_bash_file(self.test_files_path,'temp.sh'),self.js)
        os.remove(self.get_testfile_path('temp.sh'))
            
            
            
            
            
            
            
            
            
            
            