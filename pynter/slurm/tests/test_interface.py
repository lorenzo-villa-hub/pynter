#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 11:32:13 2023

@author: villa
"""
from pynter.slurm.interface import HPCInterface

from pynter.testing.core import PynterTest


class TestHPCInterface(PynterTest):
    
    def setUp(self):
        config = {'hostname': 'cluster',
                  'localdir': '/home/test/local',
                  'workdir': '/home/test/remote'}

        self.hpc = HPCInterface(config)

    def test_command(self):
        command = self.hpc.command('ls -a',printout=False,dry_run=True)[0]
        assert command == 'sshpass ssh cluster "ls" "-a" '
        
    def test_cancel_jobs(self):
        cancel_jobs = self.hpc.cancel_jobs('id-1','id-2','id-3',printout=False,dry_run=True)[0]
        assert cancel_jobs == 'sshpass ssh cluster "scancel" "id-1" "id-2" "id-3" '
        
    def test_mkdir(self):
        mkdir = self.hpc.mkdir('test',printout=False,dry_run=True)[0]
        assert mkdir == 'sshpass ssh cluster "mkdir" "-p" "test" '
        
    def test_qstat(self):
        qstat = self.hpc.qstat(printout=False,dry_run=True)[0]
        assert qstat == 'sshpass ssh cluster "squeue" "-o" ""%.10i" "%.9P" "%.40j" "%.8u" "%.2t" "%.10M" "%.5D" "%R"" '
        
    def test_sbatch(self):
        sbatch = self.hpc.sbatch('test',printout=False,dry_run=True)[0]
        assert sbatch == 'sshpass ssh cluster cd /home/test/remote/test ; sbatch job.sh'
        
        #rsync cannot be tested without connecting to a real cluster
    
    
    