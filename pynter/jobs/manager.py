
import os.path as op
from pynter import run_local
    


def get_cancel_jobs_command(*args):
    """
    Get command to cancel jobs using scancel

    Parameters
    ----------
    *args : (list)
        Job-IDs.
    
    Returns
    -------
    cmd : (str)
        Command.
    """
    cmd = 'scancel '
    for arg in args:
        cmd += arg + ' '
    return cmd            


def get_mkdir_command(path):
    """
    Get command to make a new directory if doesn't exist.

    Parameters
    ----------
    path : (str)

    Returns
    -------
    cmd : (str)
        Command.
    """
    path = op.abspath(path)
    cmd = 'mkdir -p %s' %path
    return cmd


def get_qstat_command():
    """
    Get command to check job queue status
    """
    cmd='squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R"'
    return cmd


def get_sbatch_command(self,path,job_script_filename='job.sh'):
    """
    Execute "sbatch" command to run job.

    Parameters
    ----------
    path : (str), optional
        Path where the job script is stored (relative or absolute).
    job_script_filename : (str), optional
        Filename of the job script. The default is 'job.sh'.

    Returns
    -------
    cmd (str)
        Command.
    """   
    path = op.abspath(path)
    cmd = f'cd {path} ; sbatch {job_script_filename}'
    return cmd
            


class JobManager:
    
    
    def __init__(self,job):
        """
        Class to run and manage job both remotely and/or locally.
        """
 

    def cancel_jobs(self,*args,printout=True,dry_run=False,**kwargs):
        """
        Cancel jobs using scancel

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
                
        stdout,stderr = run_local(cmd,printout,dry_run,**kwargs)
        return stdout,stderr
 
    
    def mkdir(self,path,printout=True,dry_run=False,**kwargs):
        """
        Make new directory if doesn't exist.

        Parameters
        ----------
        path : (str)
            Path (relative or absolute.
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
        path = op.abspath(path)
        cmd = 'mkdir -p %s' %path
        stdout,stderr = run_local(cmd,printout,dry_run,**kwargs)
        return stdout, stderr
        
        
    def qstat(self,cmd='squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R"',printout=True,dry_run=False,**kwargs):
        """
        Check queue status.

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
        Execute "sbatch" command to run job.

        Parameters
        ----------
        path : (str), optional
            Path where the job script is stored (relative or absolute).
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
        path = op.abspath(path)
        command = f'cd {path} ; sbatch {job_script_filename}'
        stdout,stderr = run_local(command,printout,dry_run,**kwargs) # I've used the run local because there was a probelm with the multiple command given with ;
        return stdout,stderr               # to check again if possible
        
        
        
        