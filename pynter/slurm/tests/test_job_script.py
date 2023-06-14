#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 13:46:42 2023

@author: villa
"""

import os

from pynter.slurm.job_script import ScriptHandler

from pynter.testing.core import PynterTest
from pynter.testing.slurm import JobSettingsTest


class TestScriptHandler(PynterTest):
    
    def setUp(self):

        self.script_string = (
            "#!/bin/sh\n"
            "#SBATCH -A projecttest0000\n"
            "#SBATCH --job-name=test\n"
            "#SBATCH --array=1-2%1\n"
            "#SBATCH --mail-user=test@pynter.com\n"
            "#SBATCH --mail-type=ALL\n"
            "#SBATCH --nodes=4\n"
            "#SBATCH --ntasks-per-node=96\n"
            "#SBATCH --cpus-per-task=1\n"
            "#SBATCH --output=out.%j\n"
            "#SBATCH --error=err.%j\n"
            "#SBATCH --time=24:00:00\n"
            "#SBATCH --exclusive\n"
            "#SBATCH --mem-per-cpu=3500\n"
            "\n"
            "module purge\n"
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
        
        
        self.common_settings = {
            'project_id':'projecttest0000',
            'name':'test',
            'array_size':2,
            'email':'test@pynter.com',
            'path_exe':'/home/test/code',
            'add_stop_array':True,
            'add_automation':'automation.py',
            'filename':'job_test.sh'
            }
        
        self.settings = {
            'add_automation': 'automation.py',
            'add_lines_body': ['test_BODY', 'test_BODY2'],
            'add_lines_header': ['test_HEADER', 'test_HEADER2'],
            'add_stop_array': True,
            'array_size': 2,
            'cores_per_node': 96,
            'email': 'test@pynter.com',
            'error_filename': 'err.%j',
            'filename': 'job_test.sh',
            'memory_per_cpu': 3500,
            'modules': ['intel/2020.4', 'intelmpi/2020.4', 'fftw/3.3.10'],
            'name': 'test',
            'nodes': 4,
            'output_filename': 'out.%j',
            'partition': None,
            'path_exe': '/home/test/code',
            'processor': None,
            'project_id': 'projecttest0000',
            'timelimit': '24:00:00'
             }
        
        self.sh = ScriptHandler(**self.settings)


    def test_string(self):
        assert self.sh.__str__() == self.script_string
        assert self.sh.settings == self.settings
        
    def test_from_file(self):
        sh_from_file = ScriptHandler.from_file(self.test_files_path,'job_test.sh')
        JobSettingsTest().assert_job_settings_equal(self.sh.settings,sh_from_file.settings)
        
    def test_write_script(self):
        self.sh.filename = 'temp.sh'
        self.sh.write_script(self.test_files_path)
        JobSettingsTest().assert_job_settings_equal(
            self.sh.settings,ScriptHandler.from_file(self.test_files_path,'temp.sh').settings)
        os.remove(self.get_testfile_path('temp.sh'))
            
            
            
            
            
            
            
            
            
            
            