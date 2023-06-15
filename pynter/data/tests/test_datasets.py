#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:25:45 2023

@author: villa
"""

from pynter.data.datasets import Dataset

from pynter.testing.core import PynterTest
from pynter.testing.data import DatasetTest


class TestDataset(PynterTest):
    
    def setUp(self):
        self.ds = Dataset.from_json(self.get_testfile_path('ds_Si_vacancies_pbe_relaxation.json'))

    def test_from_json(self):
        ds = self.ds
        assert ds.name == 'Vac_Si'
        assert ds.groups == ['q-1','q0','q1']
        assert ds.jobs_table().__class__.__name__ == 'DataFrame'
        
    def test_select_jobs(self):
        ds = self.ds
        DatasetTest().assert_jobs_equal(ds.select_jobs(groups=['q1']),[ds[4],ds[5]])
        DatasetTest().assert_jobs_equal(ds.select_jobs(groups=['q1'],exclude=True),[ds[i] for i in range(0,4)])
        DatasetTest().assert_jobs_equal(ds.select_jobs(groups=['q1'],common_node='1-PBE-SCF'),[ds[4]])
        DatasetTest().assert_jobs_equal(ds.select_jobs(groups=['q1'],common_node='1-PBE-SCF',mode='or'),[ds[i] for i in (0,2,4,5)])
        j1,j2 = ds[4].copy(),ds[5].copy()
        ds_test = Dataset([j1,j2])
        ds_filtered = ds.filter_jobs(groups=['q1'])
        DatasetTest().assert_dataset_equal(ds_filtered,ds_test)
        
    def test_sort_jobs(self):
        ds = self.ds
        DatasetTest().assert_jobs_equal( ds.sort_jobs(features='charge',reset=False), ds.sort_jobs(features='name',reset=False))
        sorted_jobs = [ds[i] for i in [0,2,4,1,3,5]]
        sorted_ds_test = Dataset(sorted_jobs,path=ds.path,name=ds.name,sort=False)
        ds.sort_jobs(features='nodes',reset=True)
        DatasetTest().assert_jobs_equal(ds.jobs,sorted_ds_test.jobs)
    