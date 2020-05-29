#!/usr/bin/env python

import os
import os.path as op
import argparse as ap
from pynter.__init__ import load_config
from pynter.tools.grep import grep

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
            processor: (str) CPU architecture
            modules: (list) List of modules to be loaded
            path_exe: (str) Path to executable \n             
            add_stop_array : (Bool), Add lines for stopping array jobs when calculation is converged. \n                
                If True is effective only if key 'array_size' is not None
            add_automation : (str) , Automation script to add to the file.                
            add_lines_header : (List) , Lines to add in the header part of the file.
            add_lines_body : (List) , Lines to add in the body part of the file.
        """
        
        config = load_config()
        default_settings = config['job_settings']
        
        for key,value in default_settings.items():
            setattr(self,key,value)
        
        for key, value in kwargs.items():
            if key in default_settings:
                setattr(self,key,value)
            else:
                raise Exception('"%s" is not a possible argument \nPossible arguments are: %s' %(key, self.settings.keys()))
  
              
    @property
    def settings(self):
        return self.__dict__
    
    @staticmethod
    def from_file(path,filename='job.sh'):
        
        d = {}
        file = op.join(path,filename)
        
        string = '#SBATCH -A project'
        line = grep(string,file)
        if line:
            line = line[-1]
            d['project_id'] = line.replace(string,'').replace('\n','')
        
        string = '#SBATCH --job-name='
        line = grep(string,file)
        if line:
            line = line[-1]
            d['name'] = line.replace(string,'').replace('\n','')
        
        string = '#SBATCH --array=1-'
        line = grep(string,file)
        if line:
            line = line[-1]
            d['array_size'] = int(line.replace(string,'').replace('%%1\n',''))
        
        string = '#SBATCH --mail-user='
        line = grep(string,file)
        if line:
            line = line[-1]
            d['email'] = line.replace(string,'').replace('\n','')
        
        string = '#SBATCH --nodes='
        line = grep(string,file)
        if line:
            line = line[-1]
            d['nodes'] = int(line.replace(string,'').replace('\n',''))
        
        string = '#SBATCH --ntasks-per-node='
        line = grep(string,file)
        if line:
            line = line[-1]
            d['cores_per_node'] = int(line.replace(string,'').replace('\n',''))
        
        string = '#SBATCH --output='
        line = grep(string,file)
        if line:
            line = line[-1]
            d['output_filename'] = line.replace(string,'').replace('\n','')
        
        string = '#SBATCH --error='
        line = grep(string,file)
        if line:
            line = line[-1]
            d['error_filename'] = line.replace(string,'').replace('\n','')
        
        string = '#SBATCH --time='
        line = grep(string,file)
        if line:
            line = line[-1]
            d['timelimit'] = line.replace(string,'').replace('\n','')
        
        string = '#SBATCH --mem-per-cpu='
        line = grep(string,file)
        if line:
            line = line[-1]
            d['memory_per_cpu'] = int(line.replace(string,'').replace('\n',''))
            
        string = '#SBATCH -C '
        line = grep(string,file)
        if line:
            line = line[-1]
            d['processor'] = line.replace(string,'').replace('\n','')
            
        string = 'ml '
        lines = grep(string,file)
        if lines:
            d['modules'] = []
            for line in lines:
                if list(line[-1])[0] != '#':
                    d['modules'].append(line.replace(string,'').replace(' \n',''))
                
        string = 'srun '
        line = grep(string,file)
        if line:
            line = line[-1]
            if list(line)[0] != '#':
                d['path_exe'] = line.replace(string,'').replace('\n','')
        
        string = "if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then"
        line = grep(string,file)
        if line:
            if list(line)[0] != '#':
                d['add_stop_array'] = True
            
        string = 'automation'
        line = grep(string,file)
        if line:
            line = line[-1]
            if list(line)[0] != '#':
                d['add_automation'] = line.replace('\n','')
            else:
                d['add_automation'] = None
            
        d['filename'] = filename
        
        return ScriptHandler(**d)
        
    
    
    # write here staticmethod from_file
    def args(self):
        """
        Add and parse arguments from command line
        """
        
        parser = ap.ArgumentParser()
        
        parser.add_argument('-A','--project',help='Project ID, default is 01136',required=False,default=self.project_id,type=str,metavar='',dest='project_id')
        parser.add_argument('-n','--name',help='Job name',required=False,default=self.name,type=str,metavar='',dest='name')
        parser.add_argument('-a','--array',help='Size of job array',required=False,default=self.array_size,type=int,metavar='',dest='array_size')
        parser.add_argument('-e','--email',help='Email address for job notification',required=False,default=self.email,type=str,metavar='',dest='email')
        parser.add_argument('-N','--nodes',help='Number of nodes, default is 4',required=False,default=self.nodes,type=int,metavar='',dest='nodes')
        parser.add_argument('-c','--cores-per-node',help='Number of cores per node, default is 24',required=False,default=self.cores_per_node,type=int,metavar='',dest='cores_per_node')
        parser.add_argument('-out','--output',help='Output filename',required=False,default=self.name,type=str,metavar='',dest='output_filename')
        parser.add_argument('-err','--error',help='Error filename',required=False,default=self.name,type=str,metavar='',dest='error_filename')
        parser.add_argument('-t','--timelimit',help='Timelimit, default is 24:00:00',required=False,default=self.timelimit,type=str,metavar='',dest='timelimit')
        parser.add_argument('-M','--memory-per-cpu',help='Memory per cpu, default is 2400',required=False,default=self.memory_per_cpu,type=int,metavar='',dest='memory_per_cpu')
        parser.add_argument('-C','--processor',help='(avx or avx2, default is avx2)',required=False,default=self.processor,type=str,metavar='',dest='processor')
        parser.add_argument('-ml','--modules',action='append',help="Modules to load, default are 'intel/2019.2','intel/2019.3','intelmpi/2019.3','fftw/3.3.8'" ,required=False,default=self.modules,type=str,metavar='',dest='modules')
        parser.add_argument('-x','--exe',help='Path to executable, default is "/home/lv51dypu/vasp-5-3-3"',required=False,default=self.path_exe,type=str,metavar='',dest='path_exe')
        parser.add_argument('-s','--stop-array',action='store_true',help='Add lines to stop array jobs when converged, default is False',required=False,default=False,dest='add_stop_array')
        parser.add_argument('-S','--automation',help='Script with automation to add',required=False,default=self.add_automation,type=str,metavar='',dest='add_automation')
        parser.add_argument('-H','--header',action='append',help='Add line to header part of script',required=False,default=self.add_lines_header,type=str,metavar='',dest='add_lines_header')
        parser.add_argument('-B','--body',action='append',help='Add line to body part of script',required=False,default=self.add_lines_body,type=str,metavar='',dest='add_lines_body')
        parser.add_argument('-f','--filename',help='File name, default is "job.sh"',required=False,default=self.filename,type=str,metavar='',dest='filename')
        
        args = parser.parse_args()
        
        # update settings
        for key,value in args.__dict__.items():
            setattr(self,key,value)
                
        return
            
    
    def write_script(self,path=None):
        """
        Write job script 
        
        Parameters
        ----------
        path : (str), optional
            Path to write job script to. The default is None. If None work dir is used.
        """
        
        self.write_script_header(path=path)
        self.write_script_body(path=path)
        
        return
    

    def write_script_body(self,path=None):
        """
        Write body part of the job script (part after #SBATCH commands) 
        
        Parameters
        ----------
        path : (str), optional
            Path to write job script to. The default is None. If None work dir is used.
        """
        if path:
            if not os.path.exists(path):
                os.makedirs(path)        
        complete_path = os.path.join(path,self.filename) if path else self.filename              
        with open(complete_path,'a') as f:
            
            if self.array_size:
                f.write('\n')
                f.write('if [ -f CONTCAR ]\n')
                f.write('then\n')
                f.write('    cp CONTCAR POSCAR\n') # KEEP THE TAB !)
                f.write('fi\n')
            
            f.write('\n')
            f.write('srun %s\n' %self.path_exe)

            automation_written = False
            if self.array_size:
                if self.add_stop_array:
                    f.write('\n')
                    f.write('convergence.py > convergence.txt\n')
                    f.write("if  grep -q 'Electronic convergence: True' convergence.txt  = true  && grep -q 'Ionic convergence: True' convergence.txt  = true; then\n")
                    if self.add_automation:
                        f.write('    %s\n' %self.add_automation) # KEEP THE TAB!
                        automation_written = True
                    f.write('    scancel ${SLURM_ARRAY_JOB_ID}_*\n') #KEEP THE TAB!
                    f.write('fi\n')
            if self.add_automation and automation_written is False:
                f.write('\n')
                f.write('%s\n' %self.add_automation)

            if self.add_lines_body:
                f.writelines('\n'.join([l for l in self.add_lines_body]))

        return
                                

    def write_script_header(self,path=None):
        """
        Write header part of the job script (part with #SBATCH commands and module loads) 
        
        Parameters
        ----------
        path : (str), optional
            Path to write job script to. The default is None. If None work dir is used.
        """
        if path:
            if not os.path.exists(path):
                os.makedirs(path)
        complete_path = os.path.join(path,self.filename) if path else self.filename      
        with open(complete_path,'w') as f:
            
            f.write('#!/bin/sh\n')
            f.write('#SBATCH -A project%s\n' %self.project_id)
            f.write('#SBATCH --job-name=%s\n' %self.name)
            if self.array_size:
                f.write('#SBATCH --array=1-%i%%1\n' %self.array_size)
            f.write('#SBATCH --mail-user=%s\n' %self.email)
            f.write('#SBATCH --mail-type=ALL\n')
            f.write('#SBATCH --nodes=%i\n' %self.nodes)
            f.write('#SBATCH --ntasks-per-node=%i\n' %self.cores_per_node)
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write('#SBATCH --output=%s\n' %self.output_filename)
            f.write('#SBATCH --error=%s\n' %self.error_filename)
            f.write('#SBATCH --time=%s\n' %self.timelimit)
            f.write('#SBATCH -p deflt\n')
            f.write('#SBATCH --exclusive\n')
            f.write('#SBATCH --mem-per-cpu=%i\n' %self.memory_per_cpu)
            f.write('#SBATCH -C %s\n' %self.processor)
            f.writelines([' '.join(['ml', m , '\n']) for m in self.modules])
            if self.add_lines_header:
                f.writelines('\n'.join([l for l in self.add_lines_header]))

        return         
        
# part to execute if script is used directly        
if __name__ == '__main__':
    
    s = ScriptHandler()
    # getting arguments
    s.args()
    # write file
    s.write_script()
