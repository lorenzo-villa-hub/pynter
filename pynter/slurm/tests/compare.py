#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:37:44 2023

@author: lorenzo
"""
import unittest

class CompareJobSettings(unittest.TestCase):
    
    @property
    def common_keys(self):
        common_keys = [
            'add_automation',
            'add_stop_array',
            'array_size',
            'cores_per_node',
            'error_filename',
            'filename',
            'memory_per_cpu',
            'name',
            'nodes',
            'output_filename',
            'partition',
            'processor',
            'timelimit'
            ]
        return common_keys
    
    def compare_settings(self,settings1,settings2):
        for key in self.common_keys:
            self.assertEqual(settings1[key],settings2[key])
        