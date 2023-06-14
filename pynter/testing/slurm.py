#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:37:44 2023

@author: lorenzo
"""
from pynter.testing.core import PynterTest

class JobSettingsTest(PynterTest):
    """
    Provides methods to test job settings dictionary
    """
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
    
    def assert_job_settings_equal(self,settings1,settings2):
        """
        Check only common keys to avoid system dependent differences
        """
        for key in self.common_keys:
            self.assert_object_almost_equal(settings1[key],settings2[key])
        