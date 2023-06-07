#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:33:14 2023

@author: villa
"""

import os
import os.path as op
import time
from datetime import datetime
import subprocess
import schedule
from shutil import which

from monty.dev import requires



from pynter.__init__ import load_config, run_local


def setup_hpc(subparsers):
    
    parser_hpc = subparsers.add_parser('HPC',help='Interface with High Performance Computer')
    
    subparsers_hpc = parser_hpc.add_subparsers()
    parser_sync = subparsers_hpc.add_parser('sync',help='Periodically sync data with HPC in the background with rsync')
    parser_sync = setup_sync(parser_sync)
    return


def setup_sync(parser):
    config = load_config()
    hostname = config['HPC']['hostname']
    workdir = config['HPC']['workdir']
    localdir = config['HPC']['localdir']
    
    parser.add_argument('-dry','--dry-run',action='store_true',help='Dry run with rsync (default: %(default)s)',default=False,dest='dry_run')
    parser.add_argument('-exc','--exclude',action='append',help='Files to exclude from sync (default: %(default)s)',
                        default=['core.*','WAVECAR'],metavar='',type=str,dest='exclude')
    parser.add_argument('-host','--hostname',help='ssh alias of the HPC (default: %(default)s)',default=hostname,metavar='',type=str,dest='hostname')
    parser.add_argument('-ldir','--localdir',help='Local calculation directory (default: %(default)s)',default=localdir,metavar='',type=str,dest='localdir')
    parser.add_argument('-wdir','--workdir',help='Remote calculation directory (default: %(default)s)',default=workdir,metavar='',type=str,dest='workdir')
    parser.add_argument('-p','--periodic',action='store_true',help='Sync every hour with HPC (default: %(default)s)',default=False,dest='periodic')
    parser.set_defaults(func=run_sync)
    return



@requires(which("rsync"),
      "rsync needs to be installed, you can install it with 'sudo apt-get install rsync'.")
def run_sync(args):
    command = get_command(args)
    run_command(command,periodic=args.periodic)
    return

    
def get_command(args):
    workdir = op.join(args.workdir,'')  #ensure backslash at the end
    localdir = op.join(args.localdir,'')
    hostname = args.hostname
    
    cmd = "rsync -r -uavzh " #keep the spaces
    if args.dry_run:
        cmd += "--dry-run "
    if args.exclude:
        for s in args.exclude:
            cmd += f'--exclude={s} ' 
    cmd += f"-e ssh {hostname}:{workdir} {localdir} "   
    return cmd
    
    
def run_command(command,periodic=False):
    if periodic:
        homedir = os.getenv("HOME")
        wdir = os.path.join(homedir,'.pynter','.sync_hpc')
        if not os.path.exists(wdir):
            os.makedirs(wdir)
        os.chdir(wdir)
        command = command.split() #subprocess.run() takes a list with the arguments of the command
        
        def job(command):        		
        	proc = subprocess.run(command, capture_output=True, shell=False, text=True)
        	stdout = proc.stdout
        	stderr = proc.stderr 
        	dtn = datetime.now()
        	date = dtn.date().isoformat()
        	hour = dtn.ctime().split()[3]
        	path = os.path.join(wdir,date)
        	if not os.path.exists(path):
        		os.mkdir(path)
        	
        	outfile = os.path.join(path,'sync_'+ hour +'.out')
        	errfile = os.path.join(path,'sync_'+ hour +'.err')
        
        	with open(outfile,'w') as outf, open(errfile,'w') as errf:
        		outf.write(stdout)
        		errf.write(stderr)
        
        schedule.every(60).minutes.do(job,command)
        
        while True:
        	schedule.run_pending()
        	time.sleep(30)
        
    else:
        run_local(command,printout=True)
        return
        






