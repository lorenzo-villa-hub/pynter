#!/usr/bin/env python

import re
import os
import os.path as op

from pynter import SETTINGS
from pynter.tools.utils import grep_list

config = {
    'HPC': 
          {'hostname': None,
          'localdir': None,
          'workdir': None},
    'API_KEY': None, 
    'job_settings': 
        {'sbatch_kwargs':{
                        'error':'err.%j',
                        'mail-user':None,
                        'mem-per-cpu':3500,
                        'ntasks':None,                    
                        'job-name':'no_name',
                        'output':'out.%j',
                        'partition':None,
                        'processor':None,
                        'account':'',
                        'time':'01:00:00'},
        'add_automation':None,
        'add_lines_body':None,
        'add_lines_header':None,
        'add_stop_array':False,
        'array_size':None,
        'filename':'job.sh',
        'modules':None,
        'path_exe':''}
    }

class SbatchScript:
    
    # switch to key value format for sbatch args
    def __init__(self,sbatch_kwargs={},filename='job.sh',array_size=None,modules=None,path_exe=None,
                 add_stop_array=False,add_automation=False,add_lines_header=None,
                 add_lines_body=None):
        """
        Parameters
        ----------
        **kwargs : 
            array_size: (int) Number of jobs for array \n
            modules: (list) List of modules to be loaded
            path_exe: (str) Path to executable \n             
            add_stop_array : (Bool), Add lines for stopping array jobs when calculation is converged. \n                
                If True is effective only if key 'array_size' is not None
            add_automation : (str) , Automation script to add to the file.                
            add_lines_header : (List) , Lines to add in the header part of the file.
            add_lines_body : (List) , Lines to add in the body part of the file.
        """
        
        default_settings = config['job_settings']
       # default_settings = SETTINGS['job_settings']
        
        for key,value in default_settings.items():
            setattr(self,key,value)
        

    def __str__(self):
        lines = self.script_header() + self.script_body()
        string = ''.join(lines)
        return string
    
    def __repr__(self):
        return self.__str__()
    
    def __len__(self):
        return len(self.settings)

    def __iter__(self):
        return self.settings.keys().__iter__()
    
    def __getitem__(self,key):
        return self.settings[key]
    
    def __setitem__(self,key,value):
        setattr(self,key,value)
        return
    
    def __eq__(self,other):
        if isinstance(other,str):
            return self.__str__() == other
        elif isinstance(other,ScriptHandler):
            return self.settings == other.settings
        elif isinstance(other,dict):
            return self.settings == other
        
    
    @property
    def settings(self):
        return self.__dict__

    @staticmethod
    def from_file(path,filename='job.sh'):
        """
        Create ScriptHandler object from file. cannot read added lines in header and body
        """
        d = {'filename':filename}
        file = op.join(path,filename)
        with open(file) as f:
            lines = [line.rstrip('\n') for line in f]
            
        
        string = '#SBATCH --array=1-'
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            line = re.sub("[^2-9]", "", line)
            d['array_size'] = int(line)
            
            
        string = 'ml '
        target_lines = grep_list(string,lines)
        if target_lines:
            d['modules'] = []
            for line in target_lines:
                if list(line[-1])[0] != '#':
                    mod_line = line.replace(string,'')
                    mod_line = mod_line.split(' ')[0] #remove space at the end
                    d['modules'].append(mod_line) 

        string = 'srun '
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            if list(line)[0] != '#':
                d['path_exe'] = line.replace(string,'')
        
        string = "if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then"
        line = grep_list(string,lines)
        if line:
            if list(line)[0] != '#':
                d['add_stop_array'] = True
            
        string = 'automation'
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            if list(line)[0] != '#':
                d['add_automation'] = line.lstrip() # exclude tab
            else:
                d['add_automation'] = None
        else:
            d['add_automation'] = None
            
        d['filename'] = filename
        
        return SbatchScript(**d)

            
    def script_body(self):
        """
        Body lines part of the job script (part after #SBATCH commands) 
        """
        f = []
        if self.array_size:
            f.append('\n')
            f.append('if [ ! -f POSCAR_initial ] ; then\n')
            f.append('    cp POSCAR POSCAR_initial\n')
            f.append('fi\n')
            f.append('if [ -f CONTCAR ]\n')
            f.append('then\n')
            f.append('    cp CONTCAR POSCAR\n') # KEEP THE TAB !)
            f.append('fi\n')
        
        f.append('\n')
        f.append('srun %s\n' %self.path_exe)

        automation_written = False
        if self.array_size:
            if self.add_stop_array:
                f.append('\n')
                f.append('pynter analysis vasprun --convergence > convergence.txt\n')
                f.append("if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then\n")
                if self.add_automation:
                    f.append('    %s\n' %self.add_automation) # KEEP THE TAB!
                    automation_written = True
                f.append('    scancel ${SLURM_ARRAY_JOB_ID}_*\n') #KEEP THE TAB!
                f.append('fi\n')
        if self.add_automation and automation_written is False:
            f.append('\n')
            f.append('%s\n' %self.add_automation)

        if self.add_lines_body:
            for l in self.add_lines_body:
                f.append(l+'\n')

        return f
                                

    def script_header(self):
        """
        Header lines part of the job script (part with #SBATCH commands and module loads) 
        """
        f = []
        f.append('#!/bin/sh\n')
        for key,value in self.sbatch_kwargs.items():
            printed_value = '' if value is True else '=%s' %value
            f.append(f'#SBATCH --{key}{printed_value} \n')
        f.append('\n')
        f.append('module purge\n')
        if self.modules:
            for m in self.modules:
                f.append(' '.join(['ml', m , '\n']))
        if self.add_lines_header:
            for l in self.add_lines_header:
                f.append(l+'\n')

        return f
    

    def write_script(self,path=None,filename=None):
        """
        Write job script 
        
        Parameters
        ----------
        path : (str), optional
            Path to write job script to. The default is None. If None work dir is used.
        filename : (str), optional
            Filename. If None self.filename is used. The default is None.
        """
        if path:
            if not os.path.exists(path):
                os.makedirs(path)
        complete_path = os.path.join(path,self.filename) if path else self.filename      
        with open(complete_path,'w') as f:
            f.write(self.__str__())        
        return
