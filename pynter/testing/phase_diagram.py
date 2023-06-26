#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:28:02 2023

@author: villa
"""

from pynter.testing.core import PynterTest


class ChempotsTest(PynterTest):
    """
    Provides methods to test Chempots objects
    """
    def assert_Chempots_equal(self,chempots1, chempots2, **kwargs):
        """
        Both Chempots or dictionary objects are allowed as input
        """
        dict1 = chempots1 if type(chempots1) == dict else chempots1.mu
        dict2 = chempots2 if type(chempots2) == dict else chempots2.mu
        self.assert_object_almost_equal(dict1, dict2, **kwargs)
        
        
class ReservoirsTest(ChempotsTest):
    """
    Provides methods to test Reservoirs objects
    """
    def assert_Reservoirs_equal(self,reservoirs1, reservoirs2, check_reference=True, **kwargs):
        """
        Reservoirs dictionary, reference chempots and are_chempots_delta are checked.
        Set check_reference to False to check only reservoirs dictionary.
        """
        keys1 = list(reservoirs1.keys())
        keys2 = list(reservoirs2.keys())
        self.assert_object_almost_equal(keys1, keys2, **kwargs)
        for key in reservoirs1:
            chempots1 = reservoirs1[key]
            chempots2 = reservoirs2[key]
            self.assert_Chempots_equal(chempots1, chempots2, **kwargs)
        if check_reference:
            self.assert_Chempots_equal(reservoirs1.mu_refs, reservoirs2.mu_refs, **kwargs)
            assert reservoirs1.are_chempots_delta == reservoirs2.are_chempots_delta