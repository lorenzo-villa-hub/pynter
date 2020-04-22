
import subprocess
import os
import os.path as op
from pynter.__init__ import load_config, run_local

class HPCInterface:
    
    def __init__(self):
        
        config = load_config()['HPC']
        
        for key,value in config.items():
            setattr(self, key, value)      


    def cancel_jobs(self,*args):
        
        cmd = 'scancel '
        for arg in args:
            cmd += arg + ' '
        
        stdout,stderr = self.command(cmd)

            
    def command(self,cmd):   
        
        cmd_split = cmd.split()
        arg = ''
        for c in cmd_split:
            arg += '"%s" ' %c
            
        command = 'sshpass ssh %s %s' %(self.hostname, arg)
        
        print(f'{self.hostname}: {cmd} \n')
        stdout,stderr = run_local(command)
        return stdout,stderr
                
    
    def qstat(self,cmd='squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R"'):
        
        stdout,stderr = self.command(cmd)
        return stdout, stderr
    

    def rsync_from_hpc(self,remotedir=None,localdir=None):
        
        remotedir = remotedir if remotedir else self.workdir
        localdir = localdir if localdir else self.localdir
        
        cmd = f"rsync -r -uavzh -e ssh {self.hostname}:{remotedir}/* {localdir} "
        
        stdout,stderr = run_local(cmd)
        return stdout,stderr


    def rsync_to_hpc(self,localdir=None,remotedir=None):
        
        localdir = localdir if localdir else self.localdir
        remotedir = remotedir if remotedir else self.workdir
           
        cmd = f"rsync -r -uavzh -e ssh  {localdir} {self.hostname}:{remotedir} "
        
        stdout,stderr = run_local(cmd)
        return stdout,stderr
        
   
    def sbatch(self,path='',job_script_filename='job.sh'):
        
        path = op.join(self.workdir,path)
        cmd = 'sbatch ' +  op.join(path,job_script_filename)
        
        stdout,stderr = self.command(cmd)
        return stdout,stderr
        
        
        
        