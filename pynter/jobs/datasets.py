
import sys
import os
import os.path as op
import warnings
import pandas as pd
import json
from monty.json import MontyDecoder

from pynter import run_command
from pynter.jobs.core import get_job_from_directory
from pynter.jobs.manager import get_qstat_command, JobManager
from pynter.hpc.ssh import SSHProtocol, rsync_from_hpc, rsync_to_hpc, get_path_relative_to_hpc
from pynter.hpc.slurm import JobSettings
from pynter.tools.utils import get_object_feature, select_objects, sort_objects
from pynter import SETTINGS, LOCAL_DIR, REMOTE_DIR


sort_default_attr = 'path'

def _check_job_script(job_script_filenames,files):
    """
    Check if job script names are in a list of files. Can be either a str or list of str
    """
    check = False
    if isinstance(job_script_filenames,str):    
        if job_script_filenames in files:
            check= True
            job_script_filename = job_script_filenames
        else:
            job_script_filename = None
    elif isinstance(job_script_filenames,list):
        job_script_filename = None
        for s in job_script_filenames:
            if s in files:
                check = True
                job_script_filename = s

    else:
        raise ValueError('job_script_filenames must be a string or a list of strings')

    return check,job_script_filename


def find_jobs(path,
              job_script_filenames=None,
              sort=sort_default_attr,
              load_outputs=True,
              verbose=False,
              hostname=None,
              localdir=None,
              remotedir=None,
              jobs_kwargs=None):
    """
    Find jobs in all folders and subfolders contained in path.
    The folder containing jobs are selected based on the presence of the file job_script_filename

    Parameters
    ----------
    path : (str)
        Parent directory.
    job_script_filenames : (str or list), optional
        Filename of job bash script. The default is 'job.sh'. Can also be a list of strings if multiple 
        file names are present. The default is 'job.sh'.
    sort : (str or list), optional
        Sort list of jobs by features. If False or None jobs are not sorted. The default is 'name'.
    load_outputs : (bool)
        Load job outputs. The default is True.
    verbose : (bool)
        Print importing log on screen.
    jobs_kwargs : (dict), optional
        Dictionary with job class name as keys and kwargs as values. Kwargs to be used when importing job 
        from directory for each job class.

    Returns
    -------
    jobs : (list)
        List of Job objects.

    """
    jobs = []
    job_script_filenames = job_script_filenames if job_script_filenames else JobSettings().filename 
    for root , dirs, files in os.walk(path):
        if files != []:
            check_job_script, job_script_filename = _check_job_script(job_script_filenames,files)
            if check_job_script:
                try:
                    job = get_job_from_directory(path=root,
                                                 job_script_filename=job_script_filename,
                                                 load_outputs=load_outputs,
                                                 hostname=hostname,
                                                 localdir=localdir,
                                                 remotedir=remotedir,
                                                 jobs_kwargs=jobs_kwargs)
                except KeyboardInterrupt:
                    print("Keyboard interrupt detected. Exiting gracefully.")
                    sys.exit(0)
                except Exception as e:
                    warnings.warn(f'Importing Job from {root} failed: {e}')
                    job = None
                if job:
                    if verbose:
                        print(f'{job.jobclass} imported from {job.path}')
                    jobs.append(job)
    if sort:
        if not isinstance(sort, list):
            features = [sort]
        else:
            features = features
        jobs = sort_objects(jobs, features=features)
                
    return jobs



class Dataset:
    
    def __init__(self,
                 jobs=None,
                 path=None,
                 name=None,
                 sort='path',
                 hostname=None,
                 localdir=None,
                 remotedir=None): 
        """
        Class to store sets of calculations

        Parameters
        ----------
        jobs : (list), optional
            List of Job objects belonging to the dataset. The default is None.        
        path : (str), optional
            Path of dataset directory. If None the work dir is used if jobs is None,
            else the commonpath between the jobs is used. The default is None.
        name : (str), optional
            Name to assign to dataset. The default is None. If None the folder name is used.
        sort : (str or list), optional
            Sort list of jobs by feature. If False or None jobs are not sorted. The default is 'path'.
        """
        if jobs:
            path = op.commonpath([j.path for j in jobs])
            self.path = op.abspath(path)
        else:
            self.path = op.abspath(path) if path else os.getcwd()
        self.name = name if name else op.basename(self.path)
        self.sort = sort
        self.jobs = jobs
        if jobs:
            self._group_jobs()
            if sort:
                self.sort_jobs(inplace=True,features=sort)

        self.hostname = hostname or SETTINGS['HPC']['hostname']
        self.localdir = localdir or LOCAL_DIR
        self.remotedir = localdir or REMOTE_DIR
        
        self.path_in_hpc, self.path_relative, self.is_path_on_hpc = get_path_relative_to_hpc(
                                                        path=self.path,
                                                        localdir=self.localdir,
                                                        remotedir=self.remotedir)

    def __str__(self):     
        return self.jobs_table().__str__()
            
    def __repr__(self):
        return self.__str__()
    
    def __iter__(self):
        return self.jobs.__iter__()
    
    def __getitem__(self,index):
        return self.jobs.__getitem__(index)


    def as_dict(self,**kwargs):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "path_relative":self.path_relative,
             "name":self.name,
             "jobs":[j.as_dict(**kwargs) for j in self.jobs],
             "sort":self.sort,
             "hostname":self.hostname,
             "localdir":self.localdir,
             "remotedir":self.remotedir}
        return d
    

    def to_json(self,path='',**kwargs):
        """
        Save Dataset object as json string or file.

        Parameters
        ----------
        path : (str), optional
            Path to the destination file. If '' the path is set to "self.path/self.name.json".
            If None a string is exported. 

        Returns
        -------
        d : (str)
            If path is not set a string is returned.

        """
        d = self.as_dict(**kwargs)
        if path == '':
            path = op.join(self.path,self.name+'.json')
        if path:
            with open(path,'w') as file:
                json.dump(d,file)
            return
        else:
            return d.__str__()   

    
    @classmethod
    def from_dict(cls,d):
        path = op.join(LOCAL_DIR, d['path_relative'])
        name = d['name']
        jobs = []
        for j in d['jobs']:
            jobs.append(MontyDecoder().process_decoded(j))
        sort = d['sort']
        hostname = d['hostname'] if 'hostname' in d.keys() else None
        localdir = d['localdir'] if 'localdir' in d.keys() else None
        remotedir = d['remotedir'] if 'remotedir' in d.keys() else None
        return cls(jobs=jobs,
                   path=path,
                   name=name,
                   sort=sort,
                   hostname=hostname,
                   localdir=localdir,
                   remotedir=remotedir)
    
    
    @staticmethod
    def from_json(path_or_string):
        """
        Build Dataset object from json file or string.

        Parameters
        ----------
        path_or_string : (str)
            If an existing path to a file is given the object is constructed reading the json file.
            Otherwise it will be read as a string.

        Returns
        -------
        Dataset object.

        """
        if op.isfile(path_or_string):
            with open(path_or_string) as file:
                d = json.load(file)
        else:
            d = json.loads(path_or_string)
        return Dataset.from_dict(d)
        
    
    
    @staticmethod
    def from_directory(path=None,
                       job_script_filenames='job.sh',
                       sort='path',
                       load_outputs=True,
                       verbose=False,
                       hostname=None,
                       localdir=None,
                       remotedir=None,
                       jobs_kwargs=None): 
        """
        Static method to build Dataset object from a directory. Jobs are selected based on where the job bash script
        is present. VaspJobs are selected based on where all input files are present (INCAR,KPOINTS,POSCAR,POTCAR).

        Parameters
        ----------
        path : (str)
            Parent directory of the dataset. If None the current wdir is used.
       job_script_filenames : (str or list)
            Filename of job bash script. The default is 'job.sh'. Can also be a list of strings if multiple 
            file names are present.
        sort : (str or list)
            Sort list of jobs by feature. If False or None jobs are not sorted. The default is 'path'.
        load_outputs : (bool)
            Load job outputs.
        verbose : (bool)
            Print importing log on screen.
        jobs_kwargs : (dict)
            Dictionary with job class name as keys and kwargs as values. Kwargs to be used when importing job 
            from directory for each job class.
        """
        path = path if path else os.getcwd()
        jobs = find_jobs(path,
                         job_script_filenames=job_script_filenames,
                         sort=sort,
                         load_outputs=load_outputs,
                         verbose=verbose,
                         jobs_kwargs=jobs_kwargs) 
        
        return  Dataset(path=path,jobs=jobs,sort=sort)
    
    
    def copy(self):
        jobs = [job.copy() for job in self.jobs]
        return Dataset(jobs=jobs,path=self.path,name=self.name,sort=self.sort)
    
    @property
    def groups(self):
        """
        Directory names of the first subdirectories in the dataset path.
        """
        if not hasattr(self, '_groups'):
            self._groups = self._get_groups()
        return self._groups
    
    @property
    def nodes_levels(self):
        """
        Nested lists of common nodes between jobs. The indexes coincide with the index 
        of the individual node points. 
        For example:
            job1.nodes = '/A/b/1'
            job2.nodes = '/A/c/2'
            job3.nodes = '/B/c/3'
        self.nodes_levels = [['A','B']
                             ['b','c']
                             ['1','2','3']]
        """
        if not hasattr(self, '_nodes_levels'):
            self._nodes_levels = self.get_nodes_levels()
        return self._nodes_levels
            

    def _get_groups(self):
        groups = []
        commonpath = op.commonpath([op.abspath(j.path) for j in self.jobs])
        grps = []
        for j in self.jobs:
            job_group = op.relpath(op.abspath(j.path),start=commonpath).split('/')[0]
            if job_group not in grps:
                grps.append(job_group)
        for g in grps:
            if list(g)[0] != '.':
                groups.append(g)
        return groups


    def _group_jobs(self):
        """
        Create "groups" and "nodes" attributes based on the common subpaths of the jobs.
        The first subdirectories of the parent dataset dir will be the groups. All of the following
        subpaths will be nodes
        """
        if len(self.jobs) == 1:
            job = self.jobs
            job[0].group = None
            job[0].nodes = None
            job[0].node_points = None
            warnings.warn('You are creating a Dataset with only one Job')
        path = op.abspath(self.path)
        groups = self.groups
        for group in groups:
            group_path = op.join(path,group)
            for job in self.jobs:
                jpath = op.abspath(job.path)
                if f'/{group}/' in jpath or op.basename(jpath) == group:
                    nodes = jpath.replace(op.commonpath([group_path,jpath]),'')
                    job.group = group
                    job.nodes = nodes
                    job.node_points = nodes.split('/')[1:] #nodes as a list
        return


    def add_jobs(self,jobs,regroup=True):
        """
        Add Jobs to Dataset.

        Parameters
        ----------
        jobs : (list)
            List of Job objects.
        regroup : (bool), optional
            If True Jobs in Datasets are regroupped and self.path is set to the commonpath. The default is True.
        """
        if not self.jobs:
            self.jobs = []
        for j in jobs:
            self.jobs.append(j)
        if regroup:
            paths = [j.path for j in self.jobs]
            paths.append(self.path)
            commonpath = op.commonpath(paths)
            self.regroup_jobs(path=commonpath)
        return

    
    def add_jobs_from_directory(self,path,job_script_filenames='job.sh',sort='name',
                                load_outputs=True,regroup=True,verbose=False,jobs_kwargs=None):
        """
        Add jobs to the Dataset searching all folders and subfolders contained in given path. 
        Jobs are selected based on where the job bash script is present. 
        VaspJobs are selected based on where all input files are present (INCAR,KPOINTS,POSCAR,POTCAR).

        Parameters
        ----------
        path : (str)
            Parent directory of the dataset.
       job_script_filenames : (str or list), optional
            Filename of job bash script. The default is 'job.sh'. Can also be a list of strings if multiple 
            file names are present. The default is 'job.sh'.
        sort : (str or list), optional
            Sort list of jobs by feature. If False or None jobs are not sorted. The default is 'name'.
        regroup : (bool), optional
            Regroup jobs after adding new jobs list, self.path is set to the commonpath. The default is True.
        verbose : (bool)
            Print importing log on screen.
        jobs_kwargs : (dict)
            Dictionary with job class name as keys and kwargs as values. Kwargs to be used when importing job 
            from directory for each job class.
        """        
        path = op.abspath(path)
        jobs = find_jobs(path,
                         job_script_filenames=job_script_filenames,
                         sort=sort,
                         load_outputs=load_outputs,
                         verbose=verbose,
                         jobs_kwargs=jobs_kwargs)
        self.add_jobs(jobs,False)
        if regroup:
            commonpath = op.commonpath([path,self.path])
            self.regroup_jobs(path=commonpath)
        return
 

    def create_job(self,cls,group='',nodes='',inputs=None,job_settings=None,
                   outputs=None,job_script_filename='job.sh',name=None,**kwargs):
        """
        Create Job object and add it to the dataset

        Parameters
        ----------
        cls : (str) 
            Job class.
        group : (str), optional
            Name of the group to which the job belongs. The default is ''.
        nodes : (str), optional
            Name of the nodes of the job. The default is ''.
        inputs : (dict), optional
            Dictionary with input data. The default is None.
        job_settings : (dict), optional
            Dictionary with job settings. The default is None. Documentation in SbatchScript class in slurm.job_script module
        outputs : (dict), optional
            Dictionary with output data. The default is None.
        job_script_filename : (str), optional
            Filename of job script. The default is 'job.sh'.
        name : (str)
            Name of the job. If None the name is searched in the job script.
        """
        
        path = op.join(self.path,group,nodes)        
        job = cls(path=path,inputs=inputs,job_settings=job_settings,outputs=outputs,
                       job_script_filename=job_script_filename,name=name,**kwargs)
        job.group = group
        job.nodes = nodes
        if self.jobs:
            self.jobs.append(job)
        else:
            self.jobs = [job]
        
        return
        

    def delete_jobs(self,jobs,delete_files=True,safety=True):
        
        for j in jobs:
            if delete_files:
                j.delete_job_files(safety=safety)
            self.jobs.remove(j)
            print('Job "%s" removed from Dataset'%j.name)
        return
        
            
    def filter_jobs(self,inplace=False,jobs=None,mode='and',exclude=False,names=None,groups=None,
                    common_group=None,common_node=None,function=None,**kwargs):
        """
        Function to filter jobs based on different selection criteria.
        The priority of the selection criterion follows the order of the input
        parameters. When more than one criterion is present, the arg "mode" 
        determines the selection criterion.

        Parameters
        ----------
        inplace : (bool), optional
            If True update the current object, otherwise returns a new Dataset object.
        jobs : (list), optional
            List of Jobs to search, if None the self.jobs is used. The default is None.
        mode : (str), optional
            Filtering mode, possibilities are: 'and' and 'or'. The default is 'and'. 
        exclude : (bool), optional
            Exclude the jobs satisfying the criteria instead of selecting them. The default is False.
        names : (list), optional
            List of Job names. 
        groups : (list), optional
            List of groups that jobs need to belong to.
        common_group : (str), optional
            String that needs to be present in the group.
        common_node : (str), optional
            String that needs to be present in the node.
        function : (function), optional
            Specific funtion for more complex criteria. The function needs to return a bool.
        **kwargs : (dict)
            Properties that the jobs need to satisfy. Keys are referred to attributes/methods 
            present in the relative Job class. To address more than one condition relative to
            the same attribute, use lists or tuples (e.g. charge=[0,1]).

        Returns
        -------
        dataset : (Dataset)
            Dataset object with filtered jobs. If inplace is True self.jobs is updated with
            the selected jobs.
        """
        jobs = self.select_jobs(jobs=jobs,mode=mode,exclude=exclude,names=names,groups=groups,
                                common_group=common_group,common_node=common_node,
                                function=function,**kwargs)
        
        if inplace:
            self.__init__(jobs,self.path,self.name,self.sort)
            return
        else:
            jobs = [job.copy() for job in jobs]
            return Dataset(jobs,self.path,self.name,self.sort)
        
    
    def get_jobs_inputs(self,**kwargs):
        """Read inputs for all jobs from the data stored in the respective directories"""
        for j in self.jobs:
            j.get_inputs(**kwargs)
        return

    def get_jobs_outputs(self,update_only=False,**kwargs):
        """Read output for all jobs from the data stored in respective directories"""
        for j in self.jobs:
            if update_only:
                if not j.outputs:
                    j.get_outputs(**kwargs)
            else:
                j.get_outputs(**kwargs)
        return
            
    def get_jobs_output_properties(self,**kwargs):
        """
        Get Job output attributes by reading data in outputs dictionary
        """
        for j in self.jobs:
            j.get_output_properties(**kwargs)
        return


    def get_nodes_levels(self,jobs=None):
        """
        Nested lists of common nodes between jobs. The indexes coincide with the index 
        of the individual node points. If jobs is None the whole list of jobs is used.
        For example:
            job1.nodes = '/A/b/1'
            job2.nodes = '/A/c/2'
            job3.nodes = '/B/c/3'
        self.nodes_levels = [['A','B']
                             ['b','c']
                             ['1','2','3']]
        """
        nodes_levels = []
        jobs = jobs or self.jobs
        for job in jobs:
            for index, node in enumerate(job.node_points):
                if len(nodes_levels) <= index:
                    nodes_levels.append([])
                if node not in nodes_levels[index]:
                    nodes_levels[index].append(node)
        return nodes_levels

    
    def jobs_table(self,jobs=[],status=False,display=[]):
        """
        Create a pandas DataFrame object to display the jobs in the dataset.

        Parameters
        ----------
        jobs : (list), optional
            List of jobs to display in the table. The default is []. If [] the attribute "jobs" is used.
        status : (bool), optional
            Display job status in table.
        display : (list), optional
            List of kwargs with methods in the Job class. The properties referred to these will be added
            to the table. See get_object_feature for more details. The default is [].

        Returns
        -------
        df : (DataFrame object)

        """
        jobs = jobs if jobs else self.jobs                           
        table = []
        index = []
        if status:
            stdout,stderr = self.queue(stdouts=True,printout=False)
        for j in jobs:
            index.append(j.name)
            d = {}
            d['formula'] = j.formula
            d['group'] = j.group
            d['nodes'] = j.nodes
            d['is_converged'] = j.is_converged
            if status:
                d['status'] = JobManager(j).get_status_from_queue(stdout)
            for feature in display:
                if isinstance(feature,list):
                    key = feature[0]
                    for k in feature[1:]:
                        key = key + '["%s"]'%k
                else:
                    key = feature
                d[key] = get_object_feature(j,feature)
            table.append(d)
            
        df = pd.DataFrame(table,index=index)
        df.index.name = 'job_name'
        
        return df
            

    def queue(self,stdouts=False,printout=True):
        """
        Display queue from HPC. If stdouts is True returns out and err strings.
        
        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        cmd = get_qstat_command()
        if self.is_path_on_hpc:
            stdout,stderr = run_command(cmd,printout=printout)   
        else:
            ssh = SSHProtocol(SETTINGS['HPC']['hostname'])
            stdout, stderr = ssh.run_command(cmd,printout=printout)

        if stdouts:
            return stdout,stderr
        else:
            return


    def regroup_jobs(self,path=None):
        """
        Regroup jobs from a given path. If path is None the current self.path is used.
        This method is different from _group_jobs since it will reset the self.path 
        attribute to the given path.
        """
        path = path if path else op.abspath(self.path)
        self.path = path
        self._groups = self._get_groups()
        self._nodes_levels = self.get_nodes_levels()
        gjobs = {}
        groups = self.groups
        for group in groups:
            gjobs[group] = {}
            group_path = op.join(path,group)
            for job in self.jobs:
                if group in job.path:
                    jpath = op.abspath(job.path)
                    nodes = jpath.replace(op.commonpath([group_path,jpath]),'')
                    gjobs[group][nodes] = job
                    job.group = group
                    job.nodes = nodes 
                    job.node_points = nodes.split('/')[1:]
        return
            

    def select_jobs(self,jobs=None,mode='and',exclude=False,names=None,groups=None,common_group=None,
                    common_node=None,function=None,**kwargs):
        """
        Function to filter jobs based on different selection criteria.
        The priority of the selection criterion follows the order of the input
        parameters. When more than one criterion is present, the arg "mode" 
        determines the selection criterion.

        Parameters
        ----------
        jobs : (list), optional
            List of Jobs to search, if None the self.jobs is used. The default is None.
        mode : (str), optional
            Filtering mode, possibilities are: 'and' and 'or'. The default is 'and'. 
        exclude : (bool), optional
            Exclude the jobs satisfying the criteria instead of selecting them. The default is False.
        names : (list), optional
            List of Job names. 
        groups : (list), optional
            List of groups that jobs need to belong to.
        common_group : (str), optional
            String that needs to be present in the group.
        common_node : (str), optional
            String that needs to be present in the node.
        function : (function), optional
            Specific funtion for more complex criteria. The function needs to return a bool.
        **kwargs : (dict)
            Properties that the jobs need to satisfy. Keys are referred to attributes/methods 
            present in the relative Job class. To address more than one condition relative to
            the same attribute, use lists or tuples (e.g. charge=[0,1]).

        Returns
        -------
        output_jobs : (list)
            List with selected jobs.
        """
        input_jobs = jobs.copy() if jobs else self.jobs.copy() 
        functions = []
        
        if names:
            def fnames(job):
                return job.name in names
            functions.append(fnames)
        
        if groups:
            def fgroups(job):
                return job.group in groups
            functions.append(fgroups)
        
        if common_group:
            def fcommon_group(job):
                return common_group in job.group
            functions.append(fcommon_group)
        
        if common_node:
            def fcommon_node(job):
                return common_node in job.nodes
            functions.append(fcommon_node)
        
        if function:
            functions.append(function)
            
        return select_objects(objects=input_jobs,mode=mode,exclude=exclude,
                              functions=functions,**kwargs)


    def sort_jobs(self,jobs_to_sort=None,features=['name'],inplace=False,reverse=False):
        """
        Sort jobs according to the target feature. If list of jobs is not given
        and inplace is True the attribute self.jobs is rewritten with the sorted list.

        Parameters
        ----------
        jobs_to_sort : (list), optional
            List of Job objects. If None self.jobs is used. The default is None.
        features : (str or list), optional
            Features to use to sort the list (see get_object_feature).
            Both a single variable or a list can be provided. The default is 'name'.
            Can also be a function that takes the Job object as an argument.
        inplace : (bool), optional
            Reset the self.jobs attribute if jobs_to_sort is None. The default is False.
        reverse : (bool), optional
            Reverse the order in the list. The default is False.

        Returns
        -------
        sorted_jobs : (list)
            List of sorted Job objects.
        """
        if not isinstance(features,list):
            features = [features]
        jobs = jobs_to_sort if jobs_to_sort else self.jobs
        sorted_jobs = sort_objects(jobs, features, reverse)

        if not jobs_to_sort:
            if inplace:
                self.jobs = sorted_jobs
                self.sort = features
                return
            else:
                return sorted_jobs
        else:
            return sorted_jobs


    def sync_dataset_from_hpc(self,stdouts=False,exclude=None,dry_run=False):
        """
        Sync Dataset project folder from HPC to local machine

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


    def sync_dataset_to_hpc(self,stdouts=False,exclude=None,dry_run=False):
        """
        Sync Dataset project folder from local machine to HPC

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

           
    def sync_jobs(self):
        """Sync job data from HPC to local machine"""
        for j in self.jobs:
            j.sync_from_hpc()
        self.get_jobs_outputs()
        return
        
                
    def write_jobs_input(self):
        """Write jobs inputs to files"""
        for job in self.jobs:
            job.write_input()
            
        
        