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
                        'account':'',
                        'error':'err.%j',
                        'mail-user':None,
                        'mem-per-cpu':3500,
                        'ntasks':None,                    
                        'job-name':'no_name',
                        'output':'out.%j',
                        'partition':None,
                        'processor':None,
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
    def __init__(self,filename='job.sh',array_size=None,modules=None,path_exe=None,
                 add_stop_array=False,add_automation=False,add_lines_header=None,
                 add_lines_body=None,**kwargs):
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
    #    default_settings = SETTINGS['job_settings']
        
        
        self.filename = filename if filename is not None else default_settings['filename']
        self.array_size = array_size if array_size is not None else default_settings['array_size']
        self.modules = modules if modules is not None else default_settings['modules']
        self.path_exe = path_exe if path_exe is not None else default_settings['path_exe']
        self.add_stop_array = add_stop_array if add_stop_array else default_settings['add_stop_array']
        self.add_automation = add_automation if add_automation else default_settings['add_automation']
        self.add_lines_header = add_lines_header if add_lines_header is not None else default_settings['add_lines_header']
        self.add_lines_body = add_lines_body if add_lines_body is not None else default_settings['add_lines_body']
        
        sbatch_kwargs = {**default_settings['sbatch_kwargs'], **kwargs}
        self.sbatch_kwargs = dict(sorted(sbatch_kwargs.items(), key=lambda item: item[0]))
        
        self._slurm_arguments, self._slurm_arguments_legend = read_possible_slurm_arguments()
          
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
        elif isinstance(other,SbatchScript):
            return self.settings == other.settings
        elif isinstance(other,dict):
            return self.settings == other
        
    
    @property
    def settings(self):
        return self.__dict__
    
    @property
    def slurm_arguments(self):
        return self._slurm_arguments

    @property
    def slurm_arguments_legend(self):
        return self._slurm_arguments_legend


    @staticmethod
    def from_file(path,filename='job.sh'):
        """
        Create ScriptHandler object from file. cannot read added lines in header and body
        """
        file = op.join(path,filename)
        with open(file) as f:
            input_string = f.read()
        sbatch_script = SbatchScript.from_string(input_string)
        sbatch_script.filename = filename
        return sbatch_script
        
    
    @staticmethod
    def from_string(input_string):
        
        lines = input_string.split('\n')
        
        sbatch_kwargs = {}
        string = '#SBATCH'
        target_lines = grep_list(string,lines)
        for line in target_lines:            
            words = line.split(' ')
            pair =  words[1].strip('--').split('=')
            key = pair[0]
            value = pair[1] if len(pair) > 1 else ''
            if value.isdigit():
                value = int(value)
            sbatch_kwargs.update({key:value})
            
        string = '#SBATCH --array=1-'
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            line = re.sub("[^2-9]", "", line)
            array_size = int(line)
            
            
        string = 'ml '
        target_lines = grep_list(string,lines)
        if target_lines:
            modules = []
            for line in target_lines:
                if list(line[-1])[0] != '#':
                    mod_line = line.replace(string,'')
                    mod_line = mod_line.split(' ')[0] #remove space at the end
                    modules.append(mod_line) 

        string = 'srun '
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            if list(line)[0] != '#':
                path_exe = line.replace(string,'')
        
        string = "if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then"
        line = grep_list(string,lines)
        if line:
            if list(line)[0] != '#':
                add_stop_array = True
            
        string = 'automation'
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            if list(line)[0] != '#':
                add_automation = line.lstrip() # exclude tab
            else:
                add_automation = None
        else:
            add_automation = None
            
        
        return SbatchScript(array_size=array_size,modules=modules,path_exe=path_exe,
                            add_stop_array=add_stop_array,add_automation=add_automation,
                            **sbatch_kwargs)

            
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
        if self.array_size:
            f.append('#SBATCH --array=1-%i%%1\n' %self.array_size)
        for key,value in self.sbatch_kwargs.items():
            if value is not None:
                printed_value = '' if value is True else '=%s' %value
                f.append(f'#SBATCH --{key}{printed_value}\n')
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



def read_possible_slurm_arguments():
    path = os.path.abspath(__file__).strip(os.path.basename(__file__))
    filename = 'slurm_arguments.txt'
    with open(os.path.join(path,filename),'r') as file:
        lines = file.readlines()
  
    arguments = []
    arguments_legend = {}
    for line in lines:
        line = line.strip('\n')
        elements = line.split(',')
        arg = elements[0]
        arguments.append(arg)
        if len(elements) > 1:
            arg_short = elements[1]
            arguments_legend.update({arg_short:arg})

    return arguments, arguments_legend

