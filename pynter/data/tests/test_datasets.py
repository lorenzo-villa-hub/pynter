#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:25:45 2023

@author: villa
"""

import os
import os.path as op

from pynter.data.datasets import Dataset
from pynter.data.tests.compare import CompareDatasets

homedir = os.getenv("HOME")
test_files_path = op.join(homedir,'pynter/pynter/data/tests/test_files')
def get_path(filename):
    return op.join(test_files_path,filename)


def test_dataset():
    
    ds = Dataset.from_json(get_path('ds_Si_vacancies_pbe_relaxation.json'))
    assert ds.name == 'Vac_Si'
    assert ds.groups == ['q-1','q0','q1']
    assert ds.jobs_table().__class__.__name__ == 'DataFrame'
    assert ds.select_jobs(groups=['q1']) == [ds[4],ds[5]]
    assert ds.select_jobs(groups=['q1'],exclude=True) == [ds[i] for i in range(0,4)]
    assert ds.select_jobs(groups=['q1'],common_node='1-PBE-SCF') == [ds[4]]
    assert ds.select_jobs(groups=['q1'],common_node='1-PBE-SCF',mode='or') == [ds[i] for i in (0,2,4,5)]
    j1,j2 = ds[4].copy(),ds[5].copy()
    ds_test = Dataset([j1,j2])
    ds_filtered = ds.filter_jobs(groups=['q1'])
    CompareDatasets().compare_jobs(ds_filtered,ds_test)
    
    ds = Dataset.from_json(get_path('ds_Si_vacancies_pbe_relaxation.json'))
    assert ds.sort_jobs(features='charge',reset=False) == ds.sort_jobs(features='name',reset=False)
    sorted_jobs = [ds[i] for i in [0,2,4,1,3,5]]
    sorted_ds_test = Dataset(sorted_jobs,path=ds.path,name=ds.name,sort=False)
    ds.sort_jobs(features='nodes',reset=True)
    assert ds.jobs == sorted_ds_test.jobs
    