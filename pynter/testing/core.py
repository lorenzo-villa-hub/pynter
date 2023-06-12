#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 17:55:35 2023

@author: lorenzo
"""

import unittest
from numpy.testing import assert_almost_equal, assert_array_almost_equal, assert_allclose


class PynterTest(unittest.TestCase):
    
    @staticmethod
    def assert_almost_equal(actual, desired, decimal=3, **kwargs):
        return assert_almost_equal(actual, desired, decimal=decimal, **kwargs)

    @staticmethod        
    def assert_array_almost_equal(actual, desired, decimal=3, **kwargs):
        return assert_almost_equal(actual, desired, decimal=decimal, **kwargs)
    
    