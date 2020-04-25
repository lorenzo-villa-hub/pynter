
import subprocess
import os
import os.path as op
from glob import glob
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

            
    def command(self,cmd,printout=True):   
        
        cmd_split = cmd.split()
        arg = ''
        for c in cmd_split:
            arg += '"%s" ' %c
         
        command = 'sshpass ssh %s %s' %(self.hostname, arg)
        if printout:
            print(f'{self.hostname}: {cmd} \n')
        stdout,stderr = run_local(command,printout)
        return stdout,stderr
                
    
    def qstat(self,cmd='squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R"',printout=True):
        
        stdout,stderr = self.command(cmd,printout)
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
        
        list_dir = glob(op.abspath(localdir)+'/*')
        for dir in list_dir:
         #   dirname = op.basename(dir)
            cmd = f"rsync -r -uavzh -e ssh  {dir} {self.hostname}:{remotedir} "
            stdout,stderr = run_local(cmd)

        return stdout,stderr
        
   
    def sbatch(self,path='',job_script_filename='job.sh'):
        
        path = op.join(self.workdir,path)
        cmd = f'cd {path} ; sbatch {job_script_filename}'
        command = f'sshpass ssh {self.hostname} ' + cmd
        stdout,stderr = run_local(command) # I've used the run local because there was a probelm with the multiple command given with ;
        return stdout,stderr               # to check again if possible
        
        
        
        