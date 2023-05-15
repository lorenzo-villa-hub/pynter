#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 13:46:42 2023

@author: villa
"""

import os
import os.path as op

from pynter.slurm.job_script import ScriptHandler

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/phase_diagram/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

kwargs = {
    'project_id':'projecttest0000',
    'name':'job_test.sh',
    'array_size':2,
    'email':'test@pynter.com',
    'path_exe':'/home/test/code',
    'add_stop_array':True,
    'add_automation':'automation.py',
    'add_lines_header':'test_HEADER',
    'add_lines_body':'test_BODY'
    }

def test_job_script():
    sh = ScriptHandler(**kwargs)
    