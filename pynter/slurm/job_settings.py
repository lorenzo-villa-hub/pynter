#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 15:16:02 2023

@author: villa
"""
import re
import os
import os.path as op
import copy

from monty.json import MSONable

from pynter import SETTINGS
from pynter.slurm.core import Slurm
from pynter.tools.utils import grep_list

class JobSettings(dict,MSONable):
    
    
    def __init__(self,load_slurm_defaults=True,slurm=None,filename=None,array_size=None,modules=None,
                 export=None,path_exe=None,add_stop_array=False,add_automation=False,add_lines_header=None,
                 add_lines_body=None,**kwargs):
        """
        Class to handle settings for slurm job submission through a bash script. 
        Subscriptable like a dictionary.
        
        Parameters
        ----------
        slurm : (Slurm) Slurm object. Handles all argments related to #SBATCH. 
            If None the defualt values are used, which will be updated with the
            user-defined **kwargs.
        array_size: (int) Number of jobs for array \n
        modules: (list) List of modules to be loaded
        path_exe: (str) Path to executable \n             
        add_stop_array : (Bool), Add lines for stopping array jobs when calculation is converged. \n                
            If True is effective only if key 'array_size' is not None
        add_automation : (str) , Automation script to add to the file.                
        add_lines_header : (List) , Lines to add in the header part of the file.
        add_lines_body : (List) , Lines to add in the body part of the file.
        **kwargs : dict. Additional kwargs which will be automatically assigned to Slurm.
        """
        super().__init__()
        
        settings = {}
        default_settings = SETTINGS['job_settings']
        
        if slurm:
            if type(slurm) is dict:
                slurm = Slurm(load_defaults=load_slurm_defaults,**slurm)
            elif slurm.__class__.__name__ != 'Slurm':
                raise ValueError('Slurm settings must be provided either as dictionary or as Slurm object')

        elif not slurm and 'slurm' not in kwargs.keys():
            slurm = Slurm(load_defaults=load_slurm_defaults,**kwargs)
        
        settings['slurm'] = slurm
    
        settings['filename'] = filename if filename is not None else default_settings['filename']
        settings['array_size'] = array_size if array_size is not None else default_settings['array_size']
        settings['modules'] = modules if modules is not None else default_settings['modules']
        settings['export'] = export if export is not None else default_settings['export']
        settings['path_exe'] = path_exe if path_exe is not None else default_settings['path_exe']
        settings['add_stop_array'] = add_stop_array if add_stop_array else default_settings['add_stop_array']
        settings['add_automation'] = add_automation if add_automation else default_settings['add_automation']
        settings['add_lines_header'] = add_lines_header if add_lines_header is not None else default_settings['add_lines_header']
        settings['add_lines_body'] = add_lines_body if add_lines_body is not None else default_settings['add_lines_body']

        for key,value in settings.items():
            setattr(self,key,value)
            
        self.update(settings)
    
    def __getitem__(self,key):
        if key not in self.keys():
            return self.slurm[key]
        else:
            return super().__getitem__(key)
    
    def __setitem__(self,key,value):
        args = ['slurm','filename','array_size','modules','export','path_exe','add_stop_array',
                'add_automation','add_lines_header','add_lines_body']
        if key not in args:
            self.slurm.__setitem__(key, value)
        else:
            super().__setitem__(key,value)
        setattr(self, key, value)
        return
           
    
    def as_dict(self):
        d = dict(self)
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
        return d       
    
    def copy(self):        
        return copy.deepcopy(self)
    
    @classmethod
    def from_dict(cls,d):
        return cls(**{k: v for k, v in d.items() if k not in ["@module", "@class"]})
    
    
    @staticmethod
    def from_bash_file(path,filename='job.sh'):
        """
        Create JobSettings object from file. Cannot read added lines in header and body
        """
        file = op.join(path,filename)
        with open(file) as f:
            input_string = f.read()
        job_settings = JobSettings.from_bash_script(input_string)
        job_settings['filename'] = filename
        return job_settings
        
    
    @staticmethod
    def from_bash_script(input_string):
        """
        Initialize object from a bash script
        """
        lines = input_string.split('\n')
        
        slurm = Slurm.from_bash_script(input_string)
        
        string = '#SBATCH --array=1-'
        line = grep_list(string,lines)
        array_size = None
        if line:
            line = line[-1]
            line = re.sub("[^2-9]", "", line)
            array_size = int(line)
        if array_size:
            del slurm['array']

        export = None
        string = 'export '
        target_lines = grep_list(string,lines)
        if target_lines:
            export = []
            for line in target_lines:
                if list(line[-1])[0] != '#':
                    export_line = line.replace(string,'')
                    export_line = export_line.split(' ')[0] #remove space at the end
                    export.append(export_line) 

        modules = None
        string = 'ml '
        target_lines = grep_list(string,lines)
        if target_lines:
            modules = []
            for line in target_lines:
                if list(line[-1])[0] != '#':
                    mod_line = line.replace(string,'')
                    mod_line = mod_line.split(' ')[0] #remove space at the end
                    modules.append(mod_line) 
        else:
            string = 'module load '
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
        add_stop_array = False
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
            
        
        return JobSettings(load_slurm_defaults=False,slurm=slurm,array_size=array_size,
                            modules=modules,export=export,path_exe=path_exe,add_stop_array=add_stop_array,
                            add_automation=add_automation)


    def get_bash_script(self):
        lines = self.get_bash_script_header() + self.get_bash_script_body()
        string = ''.join(lines)
        return string
            
    
    def get_bash_script_body(self):
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
                                

    def get_bash_script_header(self):
        """
        Header lines part of the job script (part with #SBATCH commands and module loads) 
        """
        f = self.slurm.get_bash_script_lines()
        if self.array_size:
            f.append('#SBATCH --array=1-%i%%1\n' %self.array_size)
        f.append('\n')
        f.append('module purge\n')
        if self.export:
            for exp in self.export:
                f.append(' '.join(['export', exp , '\n']))
        if self.modules:
            for m in self.modules:
                f.append(' '.join(['ml', m , '\n']))
        if self.add_lines_header:
            for l in self.add_lines_header:
                f.append(l+'\n')

        return f
    

    def write_bash_file(self,path=None,filename=None):
        """
        Write job script file.
        
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
        filename = filename if filename else self.filename
        complete_path = os.path.join(path,self.filename) if path else filename      
        with open(complete_path,'w') as f:
            f.write(self.get_bash_script())        
        return
    
    
