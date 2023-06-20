
import os.path as op
from shutil import which

from monty.dev import requires

from pynter import run_local, SETTINGS
    

class HPCInterface:
    
    
    def __init__(self,config=None):
        """
        Class to interface commands with HPC

        Parameters
        ----------
        config : (dict), optional
            Configuration dictionary. If None is loaded from .pynter/config.yml.
        """
        if not config:
            config = SETTINGS['HPC']
        
        for key,value in config.items():
            setattr(self, key, value)      


    def cancel_jobs(self,*args,printout=True,dry_run=False,**kwargs):
        """
        Cancel jobs on HPC using scancel

        Parameters
        ----------
        *args : (str)
            Job-id.
        printout : (bool), optional
            Write output on screen. The default is True.
        dry_run : (bool), optional
            Only return the command to be run and empty error string. The default is False.
        kwargs : (dict), optional
            Kwargs for the run_local function.
        """
        cmd = 'scancel '
        for arg in args:
            cmd += arg + ' '
                
        stdout,stderr = self.command(cmd,printout,dry_run,**kwargs)
        return stdout,stderr


    @requires(which("sshpass"),
              "sshpass needs to be installed, you can install it with 'sudo apt-get install sshpass'.")   
    
    def command(self,cmd,printout=True,dry_run=False,**kwargs):   
        """
        Run command on HPC

        Parameters
        ----------
        cmd : (str)
            Command to run.
        printout : (bool), optional
            Write output on screen. The default is True.
        dry_run : (bool), optional
            Only return the command to be run and empty error string. The default is False.
        kwargs : (dict), optional
            Kwargs for the run_local function.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        
        cmd_split = cmd.split()
        arg = ''
        for c in cmd_split:
            arg += '"%s" ' %c
         
        command = 'sshpass ssh %s %s' %(self.hostname, arg)
        if printout:
            print(f'{self.hostname}: {cmd} \n')
        stdout,stderr = run_local(command,printout,dry_run,**kwargs)
        return stdout,stderr
                

    def mkdir(self,path,printout=True,dry_run=False,**kwargs):
        """
        Make new directory in HPC if doesn't exist

        Parameters
        ----------
        path : (str)
            Path of directory.
        printout : (bool), optional
            Write output on screen. The default is True.
        dry_run : (bool), optional
            Only return the command to be run and empty error string. The default is False.
        kwargs : (dict), optional
            Kwargs for the run_local function.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        cmd = 'mkdir -p %s' %path
        stdout,stderr = self.command(cmd,printout,dry_run,**kwargs)
        return stdout, stderr
        
        
    def qstat(self,cmd='squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R"',printout=True,dry_run=False,**kwargs):
        """
        Check queue status on HPC

        Parameters
        ----------
        cmd : (str), optional
            Command to run. The default is 'squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R"'.
        printout : (bool), optional
            Write output on screen. The default is True.
        dry_run : (bool), optional
            Only return the command to be run and empty error string. The default is False.
        kwargs : (dict), optional
            Kwargs for the run_local function.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        
        stdout,stderr = self.command(cmd,printout,dry_run,**kwargs)
        return stdout, stderr
    
    
    @requires(which("rsync"),
          "rsync needs to be installed, you can install it with 'sudo apt-get install rsync'.")
    
    def rsync_from_hpc(self,remotedir=None,localdir=None,exclude=None,printout=True,dry_run=False,**kwargs):
        """
        Sync folders from HPC to local machine. The command "rsync" is used. With this function all
        the folders contained in the remote dir are synched to the local dir.

        Parameters
        ----------
        remotedir : (str), optional
            Remote directory. The default is None. If None the one in config.yml file is used.
        localdir : (str), optional
            Local directory. The default is None. If None the one in config.yml file is used.
        exclude : (list), optional
            List of files to exclude in rsync. The default is None.
        dry_run : (bool), optional
            Perform dry run in rsync with --dry-run. The default is False. The dry_run in run_local is set to False.
        kwargs : (dict), optional
            Kwargs for the run_local function.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """        
        remotedir = remotedir if remotedir else self.workdir
        localdir = localdir if localdir else self.localdir
        localdir = op.abspath(localdir)
        remotedir = op.join(remotedir,'')  #ensure backslash at the end
        localdir = op.join(localdir,'')
        
        localcmd = 'mkdir -p %s' %localdir
        run_local(localcmd)
        
        cmd = "rsync -r -uavzh " #keep the spaces
        if dry_run:
            cmd += "--dry-run "
        if exclude:
            for s in exclude:
                cmd += f'--exclude={s} ' 
        cmd += f"-e ssh {self.hostname}:{remotedir} {localdir} "

        stdout,stderr = run_local(cmd,printout=printout,dry_run=False,**kwargs)
        return stdout,stderr


    @requires(which("rsync"),
          "rsync needs to be installed, you can install it with 'sudo apt-get install rsync'.")
    
    def rsync_to_hpc(self,localdir=None,remotedir=None,exclude=None,printout=True,dry_run=False,**kwargs):
        """
        Sync folders from local machine to HPC. The command "rsync" is used. With this function all
        the folders contained in the local dir are synched to the remote dir.

        Parameters
        ----------
        localdir : (str), optional
            Local directory. The default is None. If None the one in config.yml file is used.
        remotedir : (str), optional
            Remote directory. The default is None. If None the one in config.yml file is used.
        exclude : (list), optional
            List of files to exclude in rsync. The default is None.
        dry_run : (bool), optional
            Perform dry run in rsync with --dry-run. The default is False. The dry_run in run_local is set to False.
        kwargs : (dict), optional
            Kwargs for the run_local function.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """        
        localdir = localdir if localdir else self.localdir 
        remotedir = remotedir if remotedir else self.workdir
        localdir = op.abspath(localdir)
        localdir = op.join(localdir,'')  #ensure backslash at the end
        remotedir = op.join(remotedir,'')
        
        self.mkdir(remotedir,printout=False)
        
        cmd = "rsync -r -uavzh " #keep the spaces
        if dry_run:
            cmd += "--dry-run "
        if exclude:
            for s in exclude:
                cmd += f'--exclude={s} '
        cmd += f"-e ssh  {localdir} {self.hostname}:{remotedir} "
        
        stdout,stderr = run_local(cmd,printout=printout,dry_run=False,**kwargs)

        return stdout,stderr
    
    
    @requires(which("sshpass"),
              "sshpass needs to be installed, you can install it with 'sudo apt-get install sshpass'.")      
    
    def sbatch(self,path='',job_script_filename='job.sh',printout=True,dry_run=False,**kwargs):
        """
        Execute "sbatch" command on HPC to run job.

        Parameters
        ----------
        path : (str), optional
            Path where the job script is stored. The default is ''.
        job_script_filename : (str), optional
            Filename of the job script. The default is 'job.sh'.
        dry_run : (bool), optional
            Only return the command to be run and empty error string. The default is False.
        kwargs : (dict), optional
            Kwargs for the run_local function.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """    
        path = op.join(self.workdir,path)
        cmd = f'cd {path} ; sbatch {job_script_filename}'
        command = f'sshpass ssh {self.hostname} ' + cmd
        stdout,stderr = run_local(command,printout,dry_run,**kwargs) # I've used the run local because there was a probelm with the multiple command given with ;
        return stdout,stderr               # to check again if possible
        
        
        
        