

from pynter import run_local, SETTINGS
    

class JobManager:
    
    
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
        stdout,stderr = run_local(cmd,printout,dry_run,**kwargs)
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
        stdout,stderr = run_local(cmd,printout,dry_run,**kwargs)
        return stdout, stderr
    
      
    
    def sbatch(self,path,job_script_filename='job.sh',printout=True,dry_run=False,**kwargs):
        """
        Execute "sbatch" command on HPC to run job.

        Parameters
        ----------
        path : (str), optional
            Path where the job script is stored.
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
        command = f'cd {path} ; sbatch {job_script_filename}'
        stdout,stderr = run_local(command,printout,dry_run,**kwargs) # I've used the run local because there was a probelm with the multiple command given with ;
        return stdout,stderr               # to check again if possible
        
        
        
        