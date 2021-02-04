
import yaml
import os

homedir = os.getenv("HOME")
cfgfile = os.path.join(homedir,'.pynter','vasp.yml')

def load_vasp_default(cfgfile=cfgfile):
    with open(cfgfile,"r") as ymlfile:
        return yaml.load(ymlfile,Loader=yaml.FullLoader) # add Loader to not get warning