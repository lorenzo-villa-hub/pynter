#!/usr/bin/env python

import re
import os
import os.path as op

from pynter import SETTINGS
from pynter.tools.utils import grep_list

class ScriptHandler:
    
    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        **kwargs : 
            project_id: (str) Project ID \n
            name: (str) Job name \n
            array_size: (int) Number of jobs for array \n
            email: (str) Email address for info \n
            nodes: (int) Number of nodes \n
            cores_per_node: (int) Number of cores per node \n
            output_filename: (str) Filename of output file \n
            error_filename: (str) Filename of error file \n
            timelimit: (int) Time limit in hours \n
            memory_per_cpu:(int) Memory settings
            partition (str) Partition
            processor: (str) CPU architecture
            modules: (list) List of modules to be loaded
            path_exe: (str) Path to executable \n             
            add_stop_array : (Bool), Add lines for stopping array jobs when calculation is converged. \n                
                If True is effective only if key 'array_size' is not None
            add_automation : (str) , Automation script to add to the file.                
            add_lines_header : (List) , Lines to add in the header part of the file.
            add_lines_body : (List) , Lines to add in the body part of the file.
        """
        
        default_settings = SETTINGS['job_settings']
        
        for key,value in default_settings.items():
            setattr(self,key,value)
        
        for key, value in kwargs.items():
            if key in default_settings:
                setattr(self,key,value)
            else:
                raise KeyError('"%s" is not a possible argument \nPossible arguments are: %s' %(key, self.settings.keys()))
  
    
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
            
        
        string = '#SBATCH -A '
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['project_id'] = line.replace(string,'')
        
        string = '#SBATCH --job-name='
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['name'] = line.replace(string,'')
        
        string = '#SBATCH --array=1-'
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            line = re.sub("[^2-9]", "", line)
            d['array_size'] = int(line)
            
        string = '#SBATCH --mail-user='
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['email'] = line.replace(string,'')
        
        string = '#SBATCH --nodes='
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['nodes'] = int(line.replace(string,''))
        
        string = '#SBATCH --ntasks-per-node='
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['cores_per_node'] = int(line.replace(string,''))
        
        string = '#SBATCH --output='
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['output_filename'] = line.replace(string,'')
        
        string = '#SBATCH --error='
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['error_filename'] = line.replace(string,'')
        
        string = '#SBATCH --time='
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['timelimit'] = line.replace(string,'')
        
        string = '#SBATCH --mem-per-cpu='
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['memory_per_cpu'] = int(line.replace(string,''))
        
        string = '#SBATCH -p '
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['partition'] = line.replace(string,'')      
        
        string = '#SBATCH -C '
        line = grep_list(string,lines)
        if line:
            line = line[-1]
            d['processor'] = line.replace(string,'')
            
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
        
        return ScriptHandler(**d)

            
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
        f.append('#SBATCH -A %s\n' %self.project_id)
        f.append('#SBATCH --job-name=%s\n' %self.name)
        if self.array_size:
            f.append('#SBATCH --array=1-%i%%1\n' %self.array_size)
        if self.email:
            f.append('#SBATCH --mail-user=%s\n' %self.email)
            f.append('#SBATCH --mail-type=ALL\n')
        f.append('#SBATCH --nodes=%i\n' %self.nodes)
        f.append('#SBATCH --ntasks-per-node=%i\n' %self.cores_per_node)
        f.append('#SBATCH --cpus-per-task=1\n')
        f.append('#SBATCH --output=%s\n' %self.output_filename)
        f.append('#SBATCH --error=%s\n' %self.error_filename)
        f.append('#SBATCH --time=%s\n' %self.timelimit)
        if self.partition:
            f.append('#SBATCH -p %s\n' %self.partition)
        f.append('#SBATCH --exclusive\n')
        f.append('#SBATCH --mem-per-cpu=%i\n' %self.memory_per_cpu)
        if self.processor:
            f.append('#SBATCH -C %s\n' %self.processor)
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
