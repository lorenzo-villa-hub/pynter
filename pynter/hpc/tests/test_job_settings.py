#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 13:46:42 2023

@author: villa
"""

import os

from pynter.slurm.core import Slurm
from pynter.slurm.job_settings import JobSettings

from pynter.testing.core import PynterTest
from pynter.testing.slurm import JobSettingsTest


class TestJobSettings(PynterTest):
    
    def setUp(self):

        self.script_string = (
            "#!/bin/sh\n"
            "#SBATCH --account=projecttest0000\n"
            "#SBATCH --error=err.%j\n"
            "#SBATCH --job-name=test\n"
            "#SBATCH --mail-user=test@pynter.com\n"
            "#SBATCH --mem-per-cpu=3500\n"
            "#SBATCH --nodes=1\n"
            "#SBATCH --ntasks=96\n"
            "#SBATCH --output=out.%j\n"
            "#SBATCH --time=24:00:00\n"
            "#SBATCH --array=1-2%1\n"
            "\n"
            "module purge\n"
            "export APPLES \n"
            "export ORANGES \n"
            "ml intel/2020.4 \n"
            "ml intelmpi/2020.4 \n"
            "ml fftw/3.3.10 \n"
            "test_HEADER\n"
            "test_HEADER2\n"
            "\n"
            "if [ ! -f POSCAR_initial ] ; then\n"
            "    cp POSCAR POSCAR_initial\n"
            "fi\n"
            "if [ -f CONTCAR ]\n"
            "then\n"
            "    cp CONTCAR POSCAR\n"
            "fi\n"
            "\n"
            "srun /home/test/code\n"
            "\n"
            "pynter analysis vasprun --convergence > convergence.txt\n"
            "if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then\n"
            "    automation.py\n"
            "    scancel ${SLURM_ARRAY_JOB_ID}_*\n"
            "fi\n"
            "test_BODY\n"
            "test_BODY2\n"
        )
        
        self.slurm_kwargs = {
            'ntasks': 96,
            'mail-user': 'test@pynter.com',
            'error': 'err.%j',
            'job-name': 'test',
            'output': 'out.%j',
            'account': 'projecttest0000',
            'time': '24:00:00'
             }
        
        self.js = JobSettings(filename='job_test.sh',array_size=2,
                                modules=['intel/2020.4', 'intelmpi/2020.4', 'fftw/3.3.10'],
                                export=['APPLES','ORANGES'],
                                path_exe='/home/test/code',add_stop_array=True,add_automation='automation.py',
                                add_lines_header=['test_HEADER', 'test_HEADER2'],
                                add_lines_body=['test_BODY', 'test_BODY2'],
                                **self.slurm_kwargs)

        
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
            
            
            
            
            
            
            
            
            
            
            