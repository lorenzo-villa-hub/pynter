
from abc import abstractmethod
import os
import os.path as op
import shutil
import warnings

from pynter import SETTINGS, LOCAL_DIR, REMOTE_DIR
from pynter.tools.utils import get_object_feature
from pynter.hpc.slurm import JobSettings
from pynter.hpc.ssh import rsync_from_hpc, rsync_to_hpc, get_path_relative_to_hpc


from pynter.jobs.manager import JobManager


class Job:
    
    def __init__(self,
                 path=None,
                 inputs=None,
                 job_settings=None,
                 outputs=None,
                 job_script_filename=None,
                 name=None,
                 hostname=None,
                 localdir=None,
                 remotedir=None):
        """
        Class to control and organize inputs and outputs of a generic job.

        Parameters
        ----------
        path : (str), optional
            Path where job is stored. The default is None. If None the work dir is used.
        inputs : (dict), optional
            Dictionary with input data. The default is None.
        job_settings : (dict), optional
            Dictionary with job settings. The default is None. Documentation in JobSettings class in hpc.slurm module
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
        self.job_script_filename = job_script_filename if job_script_filename else SETTINGS['job_script_filename']
        
        self.hostname = hostname or SETTINGS['HPC']['hostname']
        self.localdir = localdir or LOCAL_DIR 
        self.remotedir = remotedir or REMOTE_DIR
        
        self.path_in_hpc, self.path_relative, self.is_path_on_hpc = get_path_relative_to_hpc(
                                                        path=self.path,
                                                        localdir=self.localdir,
                                                        remotedir=self.remotedir)

        if outputs:
            self.get_output_properties()
         
        if name:
            self.name = name
        elif self.job_settings:
            self.name = self.job_settings['sbatch']['job-name']
        elif op.isfile(op.join(self.path,self.job_script_filename)):
            s = JobSettings.from_bash_file(self.path,filename=self.job_script_filename)
            self.name = s.settings['sbatch']['job-name']
        else:
            self.name = 'no_name'
            
        if not self.job_settings:
            self.job_settings = {}
        self.job_settings['sbatch']['job-name'] = self.name


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
        
    
    def add_to_asedb(self,db,
                     duplicate=False,
                     converged_only=True,
                     verbose=True,
                     properties_to_write=[]):
        """
        Add Atoms object in self.outputs to ASE database

        Parameters
        ----------
        db : (ase database)
            ASE database object (SQLite3Database).
        duplicate : (bool)
            Add job to db even if an entry with same path is already present.
        converged_only : (bool)
            Add job only if job.is_converged is True.
        verbose : (bool)
            Print info on screen.
        properties_to_write: (list)
            List of job properties to add to ASE DB.
            kwargs of the form feature=get_object_feature(job,feature) will
            be added to ase.db.write. See docs of get_object_feature for details.
        """
        printout = []
        properties_to_write.append('path')
        if not converged_only or self.is_converged:
            if 'Atoms' in self.outputs.keys():
                if duplicate or all(self.path != row.path for row in db.select()):
                    try:
                        atoms = self.outputs['Atoms']
                        kwargs={}
                        for feature in properties_to_write:
                            if isinstance(feature,list):
                                key = feature[0]
                                for k in feature[1:]:
                                    key = key + '_%s'%k
                            else:
                                key = feature
                            kwargs[key] = get_object_feature(self,feature)
                        db.write(atoms=atoms,**kwargs)
                        printout.append(f'Job added to DB {self.path}')
                    except Exception as e:
                        print(e)
                        printout.append(f'Adding Atoms object to DB failed for {self}')
                else:
                    printout.append(f'Skipping job already in DB for {self.path}')
            else:
                printout.append(f'Atoms object not in outputs dict for {self.path}')
        else:
            printout.append(f'Skipping unconverged self at {self.path}')
        
        if verbose:
            for line in printout:
                print(line)
    
    
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
        stdout,stderr = rsync_from_hpc(hostname=self.hostname,
                                       remotedir=self.path_in_hpc,
                                       localdir=localdir,
                                       exclude=exclude,
                                       dry_run=dry_run)
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
        stdout,stderr = rsync_to_hpc(hostname=self.hostname,
                                     localdir=localdir,
                                     remotedir=self.path_in_hpc,
                                     exclude=exclude,
                                     dry_run=dry_run)
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
        


def get_job_from_directory(path=None,
                           job_script_filename=None,
                           load_outputs=True,
                           hostname=None,
                           localdir=None,
                           remotedir=None,
                           **kwargs):   
    """
    Get Job object from a directory. The job type is selected based on the content of the folder.
    - If ("INCAR","KPOINTS","POSCAR","POTCAR") files are present, VaspJob is initialized.
    - If ("INCAR","KPOINTS","POTCAR") files are present and "POSCAR" is not present, VaspNEBJob is initialized.
    - If ("input.in") files are present, LammpsJob is initialized. In case input, data and log files are named
      differently than "input.in", "structure.data" and "log.lammps", set it in jobs_kwargs as:
          jobs_kwargs = {'LammpsJob':{
                                      'input_filename':<input filename>,
                                      'data_filename':<data filename>,
                                      'log_filename':<log filename>
                                      }

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
    kwargs : (dict)
        Kwargs to be used when importing job from directory for each job class.
        If only one job class is imported use the kwargs as usual.
        If more than 1 job class is imported use it as:
            JobClass1=jobclass1_kwargs, JobClass2=jobclass2_kwargs

    Returns
    -------
    jobs : (list)
        List of Job objects.
    """    
    from pynter.jobs.vasp.vasp_jobs import VaspJob, VaspNEBJob
    from pynter.jobs.lammps.lammps_jobs import LammpsJob, input_default_filename
    
    path = path if path else os.getcwd()
    job_script_filename = job_script_filename if job_script_filename else SETTINGS['job_script_filename']
    
    if kwargs:
        if all([key not in ['VaspJob','VaspNEBJob','LammpsJob'] for key in kwargs.keys()]):
            warnings.warn('Job class not specified in kwargs, only do this if you are not looping over different job classes')
    
    for root,dirs,files in os.walk(path):   
        # VASP
        if all(f in files for f in ['INCAR','KPOINTS','POSCAR','POTCAR']):
            path = op.abspath(root)
            if kwargs:
                vaspjob_kwargs = kwargs['VaspJob'] if 'VaspJob' in kwargs.keys() else kwargs
            else:
                vaspjob_kwargs = {}
            job = VaspJob.from_directory(
                        path=path,
                        job_script_filename=job_script_filename,
                        load_outputs=load_outputs,
                        hostname=hostname,
                        localdir=localdir,
                        remotedir=remotedir,
                        **vaspjob_kwargs)
            return job

        # VASP + NEB
        elif all(f in files for f in ['INCAR','KPOINTS','POTCAR']) and 'POSCAR' not in files:
            path = op.abspath(root)
            job = VaspNEBJob.from_directory(
                                        path=path,
                                        job_script_filename=job_script_filename,
                                        load_outputs=load_outputs,
                                        hostname=hostname,
                                        localdir=localdir,
                                        remotedir=remotedir)
            return job

        # LAMMPS
        else:
            if kwargs:
                lammpsjob_kwargs = kwargs['LammpsJob'] if 'LammpsJob' in kwargs.keys() else kwargs
            else:
                lammpsjob_kwargs = {}
            lammps_input_filename = lammpsjob_kwargs['input_filename'] if 'input_filename' in lammpsjob_kwargs.keys() else input_default_filename
            if lammps_input_filename in files:
                job = LammpsJob.from_directory(
                                            path=path,
                                            job_script_filename=job_script_filename,
                                            load_outputs=load_outputs,
                                            hostname=hostname,
                                            localdir=localdir,
                                            remotedir=remotedir,
                                            **lammpsjob_kwargs)
                return job