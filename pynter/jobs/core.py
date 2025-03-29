
from abc import abstractmethod
import os
import os.path as op
import shutil
import warnings

from pynter import SETTINGS
from pynter.hpc.slurm import JobSettings
from pynter.hpc.ssh import rsync_from_hpc, rsync_to_hpc

from pynter.jobs.manager import JobManager


class Job:
    
    def __init__(self,path=None,inputs=None,job_settings=None,outputs=None,job_script_filename=None,name=None):
        """
        Class to control and organize inputs and outputs of a generic job.

        Parameters
        ----------
        path : (str), optional
            Path where job is stored. The default is None. If None the work dir is used.
        inputs : (dict), optional
            Dictionary with input data. The default is None.
        job_settings : (dict), optional
            Dictionary with job settings. The default is None. Documentation in JobSettings class in slurm.job_script module
        outputs : (dict), optional
            Dictionary with output data. The default is None.
        job_script_filename : (str), optional
            Filename of job script. The default is taken from the key 'filename' in the job_settings in the config file.
        name : (str)
            Name of the job. If none the name is searched in the job script.

        """
                
        self.path = op.abspath(path) if path else os.getcwd()
        self.inputs = inputs
        self.job_settings = JobSettings(**job_settings) if job_settings else JobSettings()
        self.outputs = outputs
        self.job_script_filename = job_script_filename if job_script_filename else JobSettings().filename
        
        self._localdir = SETTINGS['HPC']['localdir']
        self._remotedir = SETTINGS['HPC']['remotedir']
        
        if op.commonpath([self._remotedir,self.path]) == self._remotedir:
            self.is_path_on_hpc = True
        elif op.commonpath([self._localdir,self.path]) == self._localdir:
            self.is_path_on_hpc = False
        else:
            self.is_path_on_hpc = False
            warnings.warn('Job path is not within local or remote directory', UserWarning)
            
        self.path_relative = op.abspath(self.path).replace(self._localdir,'')
        self.path_in_hpc = self._remotedir + self.path_relative
        
        
        if outputs:
            self.get_output_properties()
        
        
        if name:
            self.name = name
        elif self.job_settings:
            self.name = self.job_settings['slurm']['job-name']
        elif op.isfile(op.join(self.path,self.job_script_filename)):
            s = JobSettings.from_bash_file(self.path,filename=self.job_script_filename)
            self.name = s.settings['slurm']['job-name']
        else:
            self.name = 'no_name'
            
        if not self.job_settings:
            self.job_settings = {}
        self.job_settings['slurm']['job-name'] = self.name


    def __str__(self):
        jobclass = self.jobclass
        if hasattr(self,'group'):
            if self.group != '':
                printout = '%s "%s" of group "%s"' %(jobclass, self.name, self.group)
            else:
                printout = '%s "%s"' %(jobclass, self.name)
        else:
            self.group = ''
            printout = '%s "%s"' %(jobclass, self.name)
        
        return printout
    
    def __repr__(self):
        return self.__str__()
        
    @property
    def jobclass(self):
        return self.__class__.__name__

    @property
    def job_id(self):
        if not self._job_id:
            job_manager = JobManager(self)
            stdout,stderr = job_manager.qstat()
            job_id = job_manager.get_job_id_from_queue(stdout)
            self._job_id = job_id
        return self._job_id
        
        
    def cancel_job(self):
        """
        Cancel job on HPC
        """      
        JobManager(self).scancel()    
        return 


    def copy(self):
        """
        Copy Job object
        """
        path=self.path
        inputs=self.inputs.copy() if self.inputs else None
        job_settings=self.job_settings.copy() if self.job_settings else None
        outputs=self.outputs.copy() if self.outputs else None
        job_script_filename=self.job_script_filename
        name=self.name
        return self.__class__(path=path,inputs=inputs,job_settings=job_settings,
                       outputs=outputs,job_script_filename=job_script_filename,
                       name=name)
    
    
    def delete_job_files(self,safety=True):
        """
        Delete Job folder (self.path)

        Parameters
        ----------
        safety : (bool), optional
            Ask confirmation to delete job. The default is True.
        """
        if safety:
            inp = input('Are you sure you want to delete Job %s ? (y/n) : ' %self.name)
            if inp in ('y','Y'):
                shutil.rmtree(self.path)
                print('Deleted Job %s'%self.name)
            else:
                print('Job %s is left unchanged'%self.name)
        else:
            shutil.rmtree(self.path)
            print('Deleted Job %s'%self.name)
        return


    @abstractmethod
    def get_inputs(self):
        pass

    @abstractmethod
    def get_outputs(self):
        pass
    
    @abstractmethod
    def get_output_properties(self):
        pass

    @abstractmethod
    def insert_in_database(self):
        pass


    def job_queue(self):
        """
        Print job queue from HPC on screen
        
        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        stdout,stderr = JobManager(self).qstat()        
        return stdout,stderr
        

    def run_job(self,write_input=True,sync=True):
        """
        Run job on HPC. Input files are automatically written and sync to HPC is performed.

        Parameters
        ----------
        write_input : (bool), optional
            Write input file stored in "inputs" dictionary. The default is True.
        sync : (bool), optional
            Sync files to HPC before running. The default is True

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        if write_input:
            self.write_input()
        if sync and not self.is_path_on_hpc:
            self.sync_to_hpc()
        stdout, stderr = JobManager(self).sbatch()       
        return stdout,stderr
    

    def sync_from_hpc(self,stdouts=False,exclude=None,dry_run=False):
        """
        Sync job data from HPC to local machine

        Parameters
        ----------
        stdouts : (bool), optional
            Return output and error strings. The default is False.
        exclude : (list), optional
            List of files to exclude in rsync. The default is None.
        dry_run : (bool), optional
            Perform dry run in rsync with --dry-run. The default is False.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.

        """
        if self.is_path_on_hpc:
            raise ValueError('Cannot sync from remote machine, only from local')
        localdir = op.abspath(self.path)
        stdout,stderr = rsync_from_hpc(remotedir=self.path_in_hpc,localdir=localdir,exclude=exclude,dry_run=dry_run)
        if stdouts:
            return stdout,stderr
        else:
            return
        
        
    def sync_to_hpc(self,stdouts=False,exclude=None,dry_run=False):
        """
        Sync job data from local machine to HPC

        Parameters
        ----------
        stdouts : (bool), optional
            Return output and error strings. The default is False.
        exclude : (list), optional
            List of files to exclude in rsync. The default is None.
        dry_run : (bool), optional
            Perform dry run in rsync with --dry-run. The default is False.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.

        """
        if self.is_path_on_hpc:
            raise ValueError('Cannot sync from remote machine, only from local')
        localdir = op.abspath(self.path)
        stdout,stderr = rsync_to_hpc(localdir=localdir,remotedir=self.path_in_hpc,exclude=exclude,dry_run=dry_run)
        if stdouts:
            return stdout,stderr
        else:
            return


    def status(self):
        """
        Get job status from HPC. Gets queue from HPC and gets Job status from the queue.

        Returns
        -------
        status : (str)
            Job status. Possible status are 'PENDING','RUNNING','NOT IN QUEUE'.
        """
        job_manager = JobManager(self)
        stdout,stderr = job_manager.qstat()
        status = job_manager.get_status_from_queue()            
        return status
 
           
    @abstractmethod
    def write_input():
        pass
        


def get_job_from_directory(path=None,job_script_filename=None,load_outputs=True,jobs_kwargs=None):   
    """
    Get Job object from a directory. The job type is selected based on the content of the folder.
    - If ("INCAR","KPOINTS","POSCAR","POTCAR") files are present, VaspJob is initialized.
    - If ("INCAR","KPOINTS","POTCAR") files are present and "POSCAR" is not present, VaspNEBJob is initialized.

    Parameters
    ----------
    path : (str)
        Directory containing job.
    job_script_filename : (str), optional
        Filename of job script. The default is taken from the key 'filename' in the job_settings in the config file.
    sort : (str or list), optional
        Sort list of jobs by features. If False or None jobs are not sorted. The default is 'name'.
    load_outputs : (bool)
        Load job outputs. The default is True.
    jobs_kwargs : (dict), optional
        Dictionary with job class name as keys and kwargs as values. Kwargs to be used when importing job 
        from directory for each job class.

    Returns
    -------
    jobs : (list)
        List of Job objects.
    """    
    from pynter.jobs.vasp.vasp_jobs import VaspJob, VaspNEBJob
    
    path = path if path else os.getcwd()
    job_script_filename = job_script_filename if job_script_filename else JobSettings().filename
        
    for root,dirs,files in os.walk(path):      
        if all(f in files for f in ['INCAR','KPOINTS','POSCAR','POTCAR']):
            path = op.abspath(root)
            if jobs_kwargs:
                kwargs = jobs_kwargs['VaspJob'] if 'VaspJob' in jobs_kwargs.keys() else {}
            else:
                kwargs = {}
            j = VaspJob.from_directory(path,job_script_filename=job_script_filename,load_outputs=load_outputs,**kwargs)
            return j
        elif all(f in files for f in ['INCAR','KPOINTS','POTCAR']) and 'POSCAR' not in files:
            path = op.abspath(root)
            j = VaspNEBJob.from_directory(path,job_script_filename=job_script_filename,load_outputs=load_outputs)
            return j