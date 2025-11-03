#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:39:37 2023

@author: villa
"""
from setuptools import setup, find_namespace_packages

with open("README.md") as file:
    long_description = file.read()

setup(
    name='pynter2',
    version='2.0.0',
    author='Lorenzo Villa',
    description='Tools for atomistic calculations, provides features for point-defect calculations',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_namespace_packages(exclude=["pynter.*.tests",
    					      "pynter.*.*.tests",
    					      "pynter.*.*.*.tests"]),
    install_requires=[
        'ase',
        'paramiko',
        'pymatgen>=2023.3.23',
        'PyYAML',
        'schedule'        
    ],
    package_data={
        "pynter":["*.yml"],
        "pynter.hpc":["*.txt"],
        "pynter.lammps":["*.yml"],
        "pynter.vasp":["*.yml"]
        },
    extra_requires={'test':'pytest'},
    entry_points={
    "console_scripts": [
        "pynter = pynter.cli.main:main"]}
)


