#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 17:55:35 2023

@author: lorenzo
"""
import inspect
import os
import os.path as op
import unittest
import numpy as np

from pymatgen.core.structure import Structure


class PynterTest(unittest.TestCase):
    

    def assert_all_close(self,actual, desired, rtol=1e-07, atol=0, err_msg="", verbose=True):
        """
        Tests if two arrays are almost equal up to some relative or absolute tolerance.
        """
        return np.testing.assert_allclose(actual, desired, rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)


    def assert_object_almost_equal(self, actual, desired, rtol=1e-07, atol=0, err_msg='', verbose=True):
        """
        Custom assertion function to check if two objects are almost equal.
        Supports dictionaries, numbers (int, float), tuples, lists, str, bool and NoneType.
        Uses numpy.testing.assert_almost_equal for numerical comparisons.
        """
        # if type(actual) != type(desired):
        #     raise AssertionError(f"Objects are of different types: {type(actual)}, {type(desired)}")
    
        if isinstance(actual, dict):
            assert len(actual) == len(desired), "Dictionaries have different lengths."
            for key in actual:
                assert key in desired, f"Key '{key}' not found in the second dictionary."
                self.assert_object_almost_equal(actual[key], desired[key], rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)
    
        elif isinstance(actual, (int, float)):
            np.testing.assert_allclose(actual, desired, rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)
    
        elif isinstance(actual, (tuple, list)):
            assert len(actual) == len(desired), "Sequences have different lengths."
            for elem1, elem2 in zip(actual, desired):
                self.assert_object_almost_equal(elem1, elem2, rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)
                
        elif isinstance(actual, str):
              self.assert_str_content_equal(actual, desired, err_msg=err_msg, verbose=verbose)
        
        elif isinstance(actual, (bool,type(None))):
            assert actual == desired
        

        else:
            raise AssertionError(f"Unsupported object type: {type(actual)}")


    def assert_str_content_equal(self, actual, desired, err_msg="", verbose=True):
        """
        Tests if two strings are equal, ignoring things like trailing spaces, etc.
        """
        lines1 = actual.split("\n")
        lines2 = desired.split("\n")
        if len(lines1) != len(lines2):
            return False
        failed = []
        for l1, l2 in zip(lines1, lines2):
            if l1.strip() != l2.strip():
                failed.append(f"{l1} != {l2}")
        return len(failed) == 0
    
    
    @property
    def incar_settings(self):
        incar_settings = {'ALGO': 'Normal',
         'AMIX': 0.2,
         'EDIFF': 1e-06,
         'EDIFFG': -0.05,
         'ENCUT': 550,
         'IBRION': 2,
         'ISIF': 2,
         'ISMEAR': 0,
         'ISPIN': 1,
         'KPAR': 4,
         'LCHARG': '.TRUE.',
         'LORBIT': 10,
         'LREAL': '.FALSE.',
         'LWAVE': '.TRUE.',
         'NELM': 200,
         'NSW': 0,
         'SIGMA': 0.05,
         'SYSTEM': 'Si',
         '#### Default PBE: system': 'Si',
         'ISYM': 2}
        return incar_settings
    
    @property
    def job_settings(self):
        job_settings = {
         'add_automation': 'automation_vasp.py --contcar --chgcar --wavecar --check-kpoints --error-check',
         'add_lines_body': None,
         'add_lines_header': None,
         'add_stop_array': True,
         'array_size': None,
         'cores_per_node': 24,
         'email': 'test@pynter.com',
         'error_filename': 'err.%j',
         'filename': 'job.sh',
         'memory_per_cpu': 2400,
         'modules': ['intel/2019.2', 'intel/2019.3', 'intelmpi/2019.3', 'fftw/3.3.8'],
         'name': 'Si-BS_PBE-el-str_3',
         'nodes': 1,
         'output_filename': 'out.%j',
         'partition': 'deflt',
         'path_exe': '/home/vasp-5-3-3',
         'processor': 'avx2',
         'project_id': 'project0000',
         'timelimit': '00:30:00'}
        return job_settings

    @property
    def structure(self):
        structure = Structure.from_dict(
        {'@module': 'pymatgen.core.structure',
         '@class': 'Structure',
         'charge': 0,
         'lattice': {'matrix': [[3.32548851, 0.0, 1.91997169],
           [1.10849617, 3.13530064, 1.91997169],
           [0.0, 0.0, 3.83994338]],
          'pbc': (True, True, True),
          'a': 3.839943374653261,
          'b': 3.839943378813096,
          'c': 3.83994338,
          'alpha': 59.99999998977525,
          'beta': 59.99999995393976,
          'gamma': 60.00000000512866,
          'volume': 40.036809671145996},
         'sites': [{'species': [{'element': 'Si', 'occu': 1}],
           'abc': [0.875, 0.875, 0.875],
           'xyz': [3.879736595, 2.74338806, 6.719900914999999],
           'label': 'Si',
           'properties': {}},
          {'species': [{'element': 'Si', 'occu': 1}],
           'abc': [0.125, 0.125, 0.125],
           'xyz': [0.554248085, 0.39191258, 0.959985845],
           'label': 'Si',
           'properties': {}}]}
            )
        return structure
    
    @property
    def test_files_path(self):
        """
        Path with test files ("./test_files")
        """
        module_file = inspect.getmodule(self.__class__).__file__
        module_dir = module_file.strip(op.basename(module_file))
        return op.join(module_dir,'test_files')
    
    def get_testfile_path(self,filename):
        return op.join(self.test_files_path,filename)
    

    
    