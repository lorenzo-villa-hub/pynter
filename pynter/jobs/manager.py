
import os.path as op
from pynter import run_command
from pynter.hpc.ssh import SSHProtocol
from pynter import SETTINGS
from pynter.tools.utils import grep_list
         

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
    Get command to view job queue.

    Returns
    -------
    cmd : (str)
        Command.
    """
    cmd = 'squeue -o "%.10i %.9P %.40j %.8u %.2t %.10M %.5D %R"'
    return cmd

def get_sbatch_command(path,job_script_filename='job.sh'):
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
            
def get_scancel_command(*args):
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



class JobManager:
    
    
    def __init__(self,job):
        """
        Class to run and manage job both remotely and/or locally.
        """
        self.hostname = SETTINGS['HPC']['hostname']
        self.job = job


    def _run_job_command(self,cmd,printout=True,dry_run=False,**kwargs):
        """
        Execute command related to Job. If the current location is on HPC cluster
        the command is run on the local OS, otherwise is run on the remote HPC with ssh.
        """
        if self.job.is_path_on_hpc:
            stdout,stderr = run_command(cmd,printout=printout,dry_run=dry_run)   
        else:
            ssh = SSHProtocol(self.hostname)
            stdout, stderr = ssh.run_command(cmd,printout=printout,dry_run=dry_run)
        return stdout, stderr


    def get_job_id_from_queue(self,stdout):
        """
        Get job ID from the queue on HPC. qstat is not explicitly called here to be able to handle
        looping without running qstat every time.
        """        
        queue = stdout.splitlines()
        job_lines = grep_list(self.job.name,queue)
        if job_lines == []:
            raise ValueError (f'Job named "{self.name}" is not currently running or pending')
        elif len(job_lines) > 1:
            raise ValueError (f'More than one job named "{self.name}" has been found in queue:\n{stdout}')
        else:
            job_line = job_lines[0].split()
            job_id = job_line[0]
        
        return job_id


    def get_status_from_queue(self,stdout):     
        """
        Get Job status from stdout of qstat command in HPC. 
        qstat is not explicitly called here to be able to handle
        looping without running qstat every time.
    
        Returns
        -------
        status : (str)
            Job status. Possible status are 'PENDING','RUNNING','NOT IN QUEUE'.
        """
        queue = stdout.splitlines()
        job_lines = grep_list(self.job.name,queue)
        if job_lines == []:
            status = 'NOT IN QUEUE'
        elif len(job_lines) > 1:
            raise ValueError (f'More than one job named "{self.name}" has been found in queue:\n{stdout}')
        else:
            job_line = job_lines[0].split()
            status = job_line[4]
            if status == 'PD':
                status = 'PENDING'
            if status == 'R':
                status = 'RUNNING'
            if status == 'CG':
                status = 'COMPLETED'
            
        return status
    

    def mkdir(self,printout=True,dry_run=False):
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

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        cmd = get_mkdir_command(self.job.path)
        stdout, stderr = self._run_job_command(cmd,printout=printout,dry_run=dry_run)
        return stdout, stderr
        
        
    def qstat(self,printout=True,dry_run=False):
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

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """        
        cmd = get_qstat_command()
        stdout, stderr = self._run_job_command(cmd,printout=printout,dry_run=dry_run)
        return stdout, stderr


    def sbatch(self,printout=True,dry_run=False):
        """
        Execute "sbatch" command to run job.

        Parameters
        ----------
        printout : (bool), optional
            Write output on screen. The default is True.
        dry_run : (bool), optional
            Only return the command to be run and empty error string. The default is False.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """   
        if self.job.is_path_on_hpc:
            path = self.job.path
        else:
            path = self.job.path_in_hpc
        cmd = get_sbatch_command(path=path,job_script_filename=self.job.job_script_filename)
        stdout, stderr = self._run_job_command(cmd,printout=printout,dry_run=dry_run)     # I've used the run local because there was a probelm with the multiple command given with ;
        return stdout,stderr                                            # to check again if possible
        
        
    def scancel(self,printout=True,dry_run=False):
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
        """
        job_id = self.job.job_id
        cmd = get_scancel_command(job_id)
        stdout, stderr = self._run_job_command(cmd,printout=printout,dry_run=dry_run)
        return stdout,stderr        
        
    
    