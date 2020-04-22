
import yaml
import os

homedir = os.getenv("HOME")
cfgfile = os.path.join(homedir,'pynter','config.yml')

def load_config(cfgfile=cfgfile):
    with open(cfgfile,"r") as ymlfile:
        return yaml.load(ymlfile,Loader=yaml.FullLoader) # add Loader to not get warning
