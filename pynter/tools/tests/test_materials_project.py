#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 10:24:50 2023

@author: villa
"""
import os 
import os.path as op
import json

from pynter.tests.__init__ import get_structure_Si
from pynter.tools.materials_project import MPDatabase
from pynter.tools.utils import save_object_as_json, get_object_from_json

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/tools/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)

def test_mp_database():
    mpd = MPDatabase('mp-149')
    structure = mpd.get_structure()
    assert structure == get_structure_Si()
    with open(get_path('Si_docs.json')) as file:
        docs_test = json.load(file)
    assert mpd.get_data() == docs_test
