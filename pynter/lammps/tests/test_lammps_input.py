#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  9 16:53:16 2025

@author: villa
"""
import os

from pynter.lammps.inputs import LammpsInput
from pynter.testing.core import PynterTest

input_string = """
units       metal
boundary    p p p
atom_style  atomic
    
read_data    structure.data

pair_style  pace/extrapolation
pair_coeff  * * /path/to/potential/yaml /path/to/potential/asi el1 el2 el3 el4

thermo 500
timestep 0.001

# test
velocity all create 100.0 42
fix 1 all nvt temp 50 300 0.1
run 10000
unfix 1

fix 1 further test

write_data out.data
"""

input_lines = [
    "units       metal",
    "boundary    p p p",
    "atom_style  atomic",
    "read_data    structure.data",
    "pair_style  pace/extrapolation",
    "pair_coeff  * * /path/to/potential/yaml /path/to/potential/asi el1 el2 el3 el4",
    "thermo 500",
    "timestep 0.001",
    "# test",
    "velocity all create 100.0 42",
    "fix 1 all nvt temp 50 300 0.1",
    "run 10000",
    "unfix 1",
    "fix 1 further test",
    "write_data out.data"
]


class TestLammpsInput(PynterTest):
    
    def setUp(self):
        self.inp = LammpsInput(input_lines)
        
        
    def test_from_string(self):
        self.assert_object_almost_equal(self.inp.get_string(),LammpsInput.from_string(input_string).get_string())
        
    def test_get_command(self):
        target_lines = self.inp.get_command('fix 1')
        desired_lines = ["fix 1 all nvt temp 50 300 0.1","unfix 1","fix 1 further test"]
        self.assert_object_almost_equal(target_lines, desired_lines)
        
    def test_set_command(self):
        self.inp.set_command('timestep','timestep 0.005')
        target_line = self.inp.get_command('timestep')[0]
        self.assertEqual(target_line,'timestep 0.005')
        
    def test_file(self):
        try:
            self.inp.write_file('test.in')
            actual = LammpsInput.from_file('test.in').lines
            desired = LammpsInput(input_lines).lines
            assert actual == desired
        finally:
            os.remove('test.in')
        
        
    