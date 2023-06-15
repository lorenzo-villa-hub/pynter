#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 10:24:50 2023

@author: villa
"""
import json

from pynter.tools.materials_project import MPDatabase

from pynter.testing.core import PynterTest
from pynter.testing.structure import StructureTest

class TestMPDatabase(PynterTest):

    def setUp(self):
        self.mpd = MPDatabase('mp-149')

    def test_get_structure(self):
        structure = self.mpd.get_structure()
        StructureTest().assert_Structure_equal( structure , self.structure )
        
    def test_get_data(self):
        with open(self.get_testfile_path('Si_docs.json')) as file:
            docs_test = json.load(file)
        self.assert_object_almost_equal( self.mpd.get_data() , docs_test )
