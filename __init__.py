
import yaml
import os
import subprocess

homedir = os.getenv("HOME")
cfgfile = os.path.join(homedir,'pynter','config_lenovo.yml')

def load_config(cfgfile=cfgfile):
    with open(cfgfile,"r") as ymlfile:
        return yaml.load(ymlfile,Loader=yaml.FullLoader) # add Loader to not get warning

def run_local(cmd):
    command = cmd.split()
    proc = subprocess.run(command, capture_output=True, shell=False, text=True)
    stdout = proc.stdout
    stderr = proc.stderr  
    print(stdout)
    if stderr:
        print(stderr)
    return stdout, stderr