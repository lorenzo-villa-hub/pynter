#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 13:46:42 2023

@author: villa
"""

import os
import os.path as op

from pynter.slurm.job_script import ScriptHandler
from pynter.tests.__init__ import get_job_settings

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/slurm/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

script_string = ("#!/bin/sh\n#SBATCH -A projecttest0000\n#SBATCH --job-name=test\n"
                 "#SBATCH --array=1-2%1\n#SBATCH --mail-user=test@pynter.com\n"
                 "#SBATCH --mail-type=ALL\n#SBATCH --nodes=4\n#SBATCH --ntasks-per-node=96\n"
                 "#SBATCH --cpus-per-task=1\n#SBATCH --output=out.%j\n#SBATCH --error=err.%j\n"
                 "#SBATCH --time=24:00:00\n#SBATCH --exclusive\n#SBATCH --mem-per-cpu=3500\n\n"
                 "module purge\nml intel/2020.4 \nml intelmpi/2020.4 \nml fftw/3.3.10 \n"
                 "test_HEADER\ntest_HEADER2\n\nif [ ! -f POSCAR_initial ] ; then\n"
                 "    cp POSCAR POSCAR_initial\nfi\nif [ -f CONTCAR ]\nthen\n    cp CONTCAR POSCAR\nfi\n\n"
                 "srun /home/test/code\n\nconvergence.py > convergence.txt\n"
                 "if  grep -q 'Electronic convergence: True' convergence.txt  = true  "
                 "&& grep -q 'Ionic convergence: True' convergence.txt  = true; then\n  "
                 "  automation.py\n    scancel ${SLURM_ARRAY_JOB_ID}_*\nfi\ntest_BODY\ntest_BODY2\n")

kwargs = {
    'project_id':'projecttest0000',
    'name':'test',
    'array_size':2,
    'email':'test@pynter.com',
    'path_exe':'/home/test/code',
    'add_stop_array':True,
    'add_automation':'automation.py',
    'add_lines_header':['test_HEADER','test_HEADER2'],
    'add_lines_body':['test_BODY','test_BODY2'],
    'filename':'job_test.sh'
    }

settings = {
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


def test_job_script():
    sh = ScriptHandler(**kwargs)
    for key in kwargs:
        assert kwargs[key] == sh.settings[key]
    sh = ScriptHandler(**settings)
    assert sh.__str__() == script_string
    assert sh.settings == settings
    sh_from_file = ScriptHandler.from_file(test_files_path,'job_test.sh')
    sh.add_lines_header = None #from_file can't read added lines in header and body
    sh['add_lines_body'] = None
    assert sh.settings == sh_from_file.settings
    sh.filename = 'temp.sh'
    sh.write_script(test_files_path)
    assert sh == ScriptHandler.from_file(test_files_path,'temp.sh')
    os.remove(get_path('temp.sh'))