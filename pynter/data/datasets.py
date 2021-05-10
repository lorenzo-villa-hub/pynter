
import os
import os.path as op
import operator
from pynter.data.jobs import VaspJob, VaspNEBJob
from pynter.slurm.interface import HPCInterface
import pandas as pd
import importlib
import json
from glob import glob


def _check_job_script(job_script_filenames,files):
    """ Check if job script names are in a list of files. Can be either a str or list of str"""
    check = False
    if isinstance(job_script_filenames,str):    
        if job_script_filenames in files:
            check= True
            job_script_filename = job_script_filenames
        else:
            job_script_filename = None
    elif isinstance(job_script_filenames,list):
        for s in job_script_filenames:
            if s in files:
                check = True
                job_script_filename = s
            else:
                job_script_filename = None
    else:
        raise ValueError('job_script_filenames must be a string or a list of strings')

    return check,job_script_filename


def find_jobs(path,job_script_filenames='job.sh',sort_by_name=True,load_outputs=True,jobs_kwargs=None):
    """
    Find jobs in all folders and subfolders contained in path.
    The folder contained jobs are selected based on the presence of the file job_script_filename

    Parameters
    ----------
    path : (str)
        Parent directory.
    job_script_filenames : (str or list), optional
        Filename of job bash script. The default is 'job.sh'. Can also be a list of strings if multiple 
        file names are present. The default is 'job.sh'.
    sort_by_name : (bool), optional
        Sort list of jobs by attribute "name". The default is True.
    jobs_kwargs : (dict), optional
        Dictionay with job class name as keys and kwargs as values. Kwargs to be used when importing job 
        from directory for each job class.

    Returns
    -------
    jobs : (list)
        List of Job objects.

    """
    jobs = []    
    for root , dirs, files in os.walk(path):
        if files != []:
            check_job_script, job_script_filename = _check_job_script(job_script_filenames,files)
            if check_job_script:
                if all(f in files for f in ['INCAR','KPOINTS','POSCAR','POTCAR']):
                    path = op.abspath(root)
                    if jobs_kwargs:
                        kwargs = jobs_kwargs['VaspJob'] if 'VaspJob' in jobs_kwargs.keys() else {}
                    else:
                        kwargs = {}
                    j = VaspJob.from_directory(path,job_script_filename=job_script_filename,load_outputs=load_outputs,**kwargs)
                    j.job_script_filename = job_script_filename
                    jobs.append(j)
                elif all(f in files for f in ['INCAR','KPOINTS','POTCAR']) and 'POSCAR' not in files:
                    path = op.abspath(root)
                    j = VaspNEBJob.from_directory(path,job_script_filename=job_script_filename,load_outputs=load_outputs)
                    j.job_script_filename = job_script_filename
                    jobs.append(j)
    if sort_by_name:
        jobs = sorted(jobs, key=operator.attrgetter('name'))
                
    return jobs



class Dataset:
    
    def __init__(self,jobs=None,path=None,name=None,sort_by_name=True): 
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
        sort_by_name : (bool), optional
            Sort list of jobs based on Job names. The default is True.
        """
        if jobs:
            path = op.commonpath([j.path for j in jobs])
            self.path = op.abspath(path)
        else:
            self.path = op.abspath(path) if path else os.getcwd()
        self.name = name if name else op.basename(self.path)
        self.sort_by_name = sort_by_name
        self.jobs = jobs
        if jobs:
            self._group_jobs()
            if sort_by_name:
                self.jobs = sorted(self.jobs, key=operator.attrgetter('name'))

        self._localdir = HPCInterface().localdir
        self._workdir = HPCInterface().workdir
        self._path_relative = self.path.replace(self._localdir,'')
        
        self.path_in_hpc = self._workdir + self._path_relative


    def __str__(self):     
        return self.jobs_table().__str__()
            
    def __repr__(self):
        return self.__str__()
    
    def __iter__(self):
        return self.jobs.__iter__()


    def as_dict(self,**kwargs):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "path":self.path,
             "name":self.name,
             "jobs":[j.as_dict(**kwargs) for j in self.jobs],
             "sort_by_name":self.sort_by_name}
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
        path = d['path']
        name = d['name']
        jobs = []
        for j in d['jobs']:
            job_class = j['@class']
            module = importlib.import_module(j['@module'])
            jobclass = getattr(module, job_class)
            jobs.append(jobclass.from_dict(j))
        sort_by_name = d['sort_by_name']
        
        return cls(jobs,path,name,sort_by_name)
    
    
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
            d = json.load(path_or_string)
        return Dataset.from_dict(d)
        
    
    
    @staticmethod
    def from_directory(path=None,job_script_filenames='job.sh',sort_by_name=True,load_outputs=True,jobs_kwargs=None): 
        """
        Static method to build Dataset object from a directory. Jobs are selected based on where the job bash script
        is present. VaspJobs are selected based on where all input files are present (INCAR,KPOINTS,POSCAR,POTCAR).

        Parameters
        ----------
        path : (str)
            Parent directory of the dataset. If None the current wdir is used.
       job_script_filenames : (str or list), optional
            Filename of job bash script. The default is 'job.sh'. Can also be a list of strings if multiple 
            file names are present. The default is 'job.sh'.
        sort_by_name : (bool), optional
            Sort list of jobs by attribute "name". The default is True.
        jobs_kwargs : (dict), optional
            Dictionay with job class name as keys and kwargs as values. Kwargs to be used when importing job 
            from directory for each job class.
        """
        path = path if path else os.getcwd()
        jobs = find_jobs(path,job_script_filenames=job_script_filenames,sort_by_name=False,
                         load_outputs=load_outputs,jobs_kwargs=jobs_kwargs) # names are sorted in __init__ method
        
        return  Dataset(path=path,jobs=jobs,sort_by_name=sort_by_name)
    
    
    @property
    def groups(self):
        """Directory names of the first subdirectories in the dataset path."""
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

    
    def add_jobs_from_directory(self,path,job_script_filenames='job.sh',sort_by_name=True,load_outputs=True,regroup=True):
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
        sort_by_name : (bool), optional
            Sort list of jobs by attribute "name". The default is True.
        regroup : (bool), optional
            Regroup jobs after adding new jobs list, self.path is set to the commonpath. The default is True.
        """        
        path = op.abspath(path)
        jobs = find_jobs(path,job_script_filenames=job_script_filenames,sort_by_name=sort_by_name,load_outputs=load_outputs)
        self.add_jobs(jobs,False)
        if regroup:
            commonpath = op.commonpath([path,self.path])
            self.regroup_jobs(path=commonpath)
        return
 

    def create_job(self,cls,group='',nodes='',inputs=None,job_settings=None,
                   outputs=None,job_script_filename='job.sh',name=None):
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
            Dictionary with job settings. The default is None. Documentation in ScriptHandler class in slurm.job_script module
        outputs : (dict), optional
            Dictionary with output data. The default is None.
        job_script_filename : (str), optional
            Filename of job script. The default is 'job.sh'.
        name : (str)
            Name of the job. If None the name is searched in the job script.
        """
        
        path = op.join(self.path,group,nodes)        
        job = cls(path=path,inputs=inputs,job_settings=job_settings,outputs=outputs,
                       job_script_filename=job_script_filename,name=name)
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
        

    def get_job_feature(self,job,feature):
        """
        Get value of attribute or method of a target Job.
        If feature is a single method only the string with the method's name is required.
        If the target feature is stored in a dictionary (or dict of dictionaries), a list of this format needs to be provided:
            ["method_name",key1,key2,...] - This will identify the value of Job.method[key1][key2][...] .

        Parameters
        ----------
        job : (Job)
            Job object.
        feature : (str or list)
            Method or attribute of Job class for which the value is needed.
        """
        if isinstance(feature,list):
            met = feature[0]
            try:
                attr = getattr(job,met) ()
            except:
                attr = getattr(job,met)                
            for k in feature[1:]:
                if isinstance(attr[k],dict):
                    attr = attr[k]
                else:
                    return attr[k]
                
        else:
            met = feature
            try:
                attr = getattr(job,met) ()
            except:
                attr = getattr(job,met)
            return attr


    def get_jobs_inputs(self):
        """Read inputs for all jobs from the data stored in the respective directories"""
        for j in self.jobs:
            j.get_inputs()
        return

    def get_jobs_outputs(self,update_only=False):
        """Read output for all jobs from the data stored in respective directories"""
        for j in self.jobs:
            if update_only:
                if not j.outputs:
                    j.get_outputs()
            else:
                j.get_outputs()
        return
            
    def get_jobs_output_properties(self,**kwargs):
        """
        Get Job output attributes by reading data in outputs dictionary
        """
        for j in self.jobs:
            j.get_output_properties(**kwargs)
        return

    
    def jobs_table(self,jobs=[],display=[]):
        """
        Create a pandas DataFrame object to display the jobs in the dataset.

        Parameters
        ----------
        jobs : (list), optional
            List of jobs to display in the table. The default is []. If [] the attribute "jobs" is used.
        display : (list), optional
            List of kwargs with methods in the Job class. The properties referred to these will be added
            to the table. See self.get_job_feature for more details. The default is [].

        Returns
        -------
        df : (DataFrame object)

        """
        jobs = jobs if jobs else self.jobs                           
        table = []
        index = []
        for j in jobs:
            index.append(j.name)
            d = {}
            d['formula'] = j.formula
            d['group'] = j.group
            d['nodes'] = j.nodes
            d['is_converged'] = j.is_converged
            for feature in display:
                if isinstance(feature,list):
                    key = feature[0]
                    for k in feature[1:]:
                        key = key + '["%s"]'%k
                else:
                    key = feature
                d[key] = self.get_job_feature(j,feature)
            table.append(d)
            
        df = pd.DataFrame(table,index=index)
        df.index.name = 'job_name'
        
        return df
            

    def insert_jobs_in_database(self,safety=True,check_convergence=True,**kwargs):
        """
        Insert all jobs into relative database based on their class.

        Parameters
        ----------
        safety : (bool), optional
            Ask confirmation to insert every job in database. The default is True.
        check_convergence: (bool), optional
            Insert jobs in DB only if is_converged is True. The default is True.
        **kwargs :
            Args to pass to specific Job Drone.
        """
        for j in self.jobs:
            j.insert_in_database(safety=safety,check_convergence=check_convergence,**kwargs)
        return


    def queue(self,stdouts=False):
        """
        Display queue from HPC. If stdouts is True returns out and err strings.
        
        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        hpc = HPCInterface()
        stdout,stderr = hpc.qstat()
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
        return


    def sort_jobs(self,jobs_to_sort=None,feature='name',reset=True,reverse=False):
        """
        Sort jobs according to the target feature. If list of jobs is not given
        and reset is True the attribute self.jobs is rewritten with the sorted list.

        Parameters
        ----------
        jobs_to_sort : (list), optional
            List of Job objects. If None self.jobs is used. The default is None.
        feature : (str or list), optional
            Feature to use to sort the list (see self.get_job_feature). The default is 'name'.
        reset : (bool), optional
            Reset the self.jobs attribute if jobs_to_sort is None. The default is True.
        reverse : (bool), optional
            Reverse the order in the list. The default is False.

        Returns
        -------
        sorted_jobs : (list)
            List of sorted Job objects.
        """
        jobs = jobs_to_sort if jobs_to_sort else self.jobs
        sorted_jobs = sorted(jobs, key=lambda x: self.get_job_feature(x,feature),reverse=reverse)

        if not jobs_to_sort:
            if reset:
                self.jobs = sorted_jobs
                return
            else:
                return sorted_jobs
        else:
            return sorted_jobs
            

    def select_jobs(self,jobs=None,names=None,groups=None,common_group=None,common_node=None,complex_features=None,**kwargs):
        """
        Function to filter jobs based on different selection criteria.
        The priority of the selection criterion follows the order of the input
        parameters. When more than one criterion is present, all of them need to be 
        satisfied for the job to be selected.

        Parameters
        ----------
        jobs : (list), optional
            List of Jobs to search, if None the self.jobs is used. The default is None.
        names : (str), optional
            Job name. The default is None.
        groups : (list), optional
            List of groups that jobs need to belong to. The default is None.
        common_group : (str), optional
            String that needs to be present in the group. The default is None.
        common_node : (str), optional
            String that needs to be present in the node. The default is None.
        complex_features : (list of tuples) , optional
            List of properties that need to be satisfied.
            To use when the property of interest is stored in a dictionary.
            The first element of the tuple indentifies the property (see self.get_job_feature),
            the second corrisponds to the target value.
        **kwargs : (dict)
            Properties that the jobs need to satisfy. Keys are referred to methods 
            present in the relative Job class.

        Returns
        -------
        sel_jobs : (list)
            List with filtered jobs. If list has one element only the element is returned.
        """
        sel_jobs = jobs.copy() if jobs else self.jobs.copy() 
        jobs = sel_jobs.copy()
        
        if names:
            for j in jobs:
                if j.name not in names:
                    sel_jobs.remove(j)
        
        if groups:
            jobs = sel_jobs.copy()
            for j in jobs:
                if j.group not in groups:
                    sel_jobs.remove(j)

        if common_group:
            jobs = sel_jobs.copy()
            for j in jobs:
                if common_group not in j.group:
                    sel_jobs.remove(j)

        if common_node:
            jobs = sel_jobs.copy()
            for j in jobs:
                if common_node not in j.nodes:
                    sel_jobs.remove(j)
        
        if complex_features:
            for feature in complex_features:
                feature_name = feature[0]
                feature_value = feature[1]
                jobs = sel_jobs.copy()
                for j in jobs:
                    job_feature = self.get_job_feature(j,feature_name)
                    if job_feature != feature_value:
                        sel_jobs.remove(j)

        for feature in kwargs:
            jobs = sel_jobs.copy()
            for j in jobs:
                job_feature = self.get_job_feature(j,feature)
                if job_feature != kwargs[feature]:
                    sel_jobs.remove(j)
         
        if len(sel_jobs) == 1:
            sel_jobs = sel_jobs[0]
            
        return sel_jobs


    def sync_dataset_from_hpc(self,stdouts=False):
        """
        Sync Dataset project folder from HPC to local machine

        Parameters
        ----------
        stdouts : (bool), optional
            Return output and error strings. The default is False.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        hpc = HPCInterface()
        abs_path = op.abspath(self.path)
        localdir = abs_path 
        stdout,stderr = hpc.rsync_from_hpc(localdir=localdir,remotedir=self.path_in_hpc)
        if stdouts:
            return stdout,stderr
        else:
            return


    def sync_dataset_to_hpc(self,stdouts=False):
        """
        Sync Dataset project folder from local machine to HPC

        Parameters
        ----------
        stdouts : (bool), optional
            Return output and error strings. The default is False.

        Returns
        -------
        stdout : (str)
            Output.
        stderr : (str)
            Error.
        """
        hpc = HPCInterface()
        abs_path = op.abspath(self.path)
        localdir = abs_path 
        stdout,stderr = hpc.rsync_to_hpc(localdir=localdir,remotedir=self.path_in_hpc)
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
            
        
        