#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:39:37 2023

@author: villa
"""

from setuptools import setup, find_namespace_packages

setup(
    name='pynter',
    version='1.0.0',
    packages=find_namespace_packages(exclude=["pynter.tests","pynter.*.tests", "pynter.*.*.tests"]),
    install_requires=[
        'ase',
        'pymatgen',
        'pymatgen-analysis-defects',
        'pymatgen-db',
        'PyYAML',
        'schedule'        
    ],
    extra_requires={'test':'pytest'}
)
