#!/usr/bin/env python

import os
import time
from glob import glob
from datetime import datetime
import subprocess
import schedule

config = {
	'hostname': 'lichtenberg',
	'workdir': '/work/scratch/lv51dypu',
	'localdir': '/nfshome/villa/local-data'
	}

hostname = config['hostname']
workdir = config['workdir']
localdir = config['localdir']

homedir = os.getenv("HOME")
wdir = os.path.join(homedir,'.sync_hpc')
if not os.path.exists(wdir):
    os.makedirs(wdir)
os.chdir(wdir)

command = f"rsync -r -uavzh -e ssh {hostname}:{workdir}/*  {localdir}"
command = command.split() #subprocess.run() takes a list with the arguments of the command

def job(command):
		
	proc = subprocess.run(command, capture_output=True, shell=False, text=True)
	stdout = proc.stdout
	stderr = proc.stderr 

	ctime = datetime.now().ctime().split()
	ctime.pop(0)
	ctime = ctime[-1:] + ctime[:-1]
	ctime = '_'.join(ctime)

	outfile = 'sync_out_'+ ctime +'.log'
	errfile = 'sync_err_'+ ctime +'.log'

	with open(outfile,'w') as outf, open(errfile,'w') as errf:
		outf.write(stdout)
		errf.write(stderr)


schedule.every(15).minutes.do(job,command)


while True:

	schedule.run_pending()
	time.sleep(10)








