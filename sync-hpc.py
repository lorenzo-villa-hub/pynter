#!/usr/bin/env python

import os
import os.path as op
import time
from datetime import datetime
import subprocess
import schedule
from pynter.__init__ import load_config

config = load_config()
hostname = config['HPC']['hostname']
workdir = config['HPC']['workdir']
localdir = config['HPC']['localdir']
workdir = op.join(workdir,'')  #ensure backslash at the end
localdir = op.join(localdir,'')

homedir = os.getenv("HOME")
wdir = os.path.join(homedir,'.sync_hpc')
if not os.path.exists(wdir):
    os.makedirs(wdir)
os.chdir(wdir)

command = f"rsync -r -uavh --exclude='core.*' --exclude='WAVECAR' -e ssh {hostname}:{workdir}  {localdir}"
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








