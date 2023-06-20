
import yaml
import os
import subprocess
import warnings

def get_config_from_default_file():
    """
    Get initial settings from config_default.yml
    """
    path = os.path.abspath(__file__).strip(os.path.basename(__file__))
    filename = 'config_default.yml'
    with open(os.path.join(path,filename),'r') as ymlfile:
        return yaml.load(ymlfile,Loader=yaml.FullLoader)
    
def get_vasp_defaults_from_default_file():
    """
    Get initial vasp settings from ./vasp/vasp_default.yml
    """
    from pynter import vasp
    path = os.path.abspath(vasp.__file__).strip(os.path.basename(vasp.__file__))
    filename = 'vasp_default.yml'
    with open(os.path.join(path,filename),'r') as ymlfile:
        return yaml.load(ymlfile,Loader=yaml.FullLoader)


SETTINGS = get_config_from_default_file()
SETTINGS['vasp'] = get_vasp_defaults_from_default_file()


def get_cfgfile():
    homedir = os.getenv("HOME")
    cfgfile = os.path.join(homedir,'.pynter','config.yml')
    return cfgfile

cfgfile = get_cfgfile()
def load_config(cfgfile=cfgfile):
    """
    Load dictionary with configuration from yaml cofiguration file. The default is in ~/.pynter/config.yml.

    Parameters
    ----------
    cfgfile : (str), optional
        Path to configuration file. The default is ~/.pynter/config.yml.

    Returns
    -------
    (dict)
        Configuration dictionary.
    """
    if os.path.exists(cfgfile):
        with open(cfgfile,"r") as ymlfile:
            return yaml.load(ymlfile,Loader=yaml.FullLoader) # add Loader to not get warning
    else:
        warnings.warn('%s does not exist. Run "pynter configure" in the terminal to create it. Using default settings for now. '%cfgfile)
        return


def get_vasp_cfgfile():
    homedir = os.getenv("HOME")
    cfgfile = os.path.join(homedir,'.pynter','vasp.yml')
    return cfgfile


cfgfile = get_vasp_cfgfile()
def load_vasp_default(cfgfile=cfgfile):
    """
    Load dictionary with VASP configuration from yaml cofiguration file. The default is in ~/.pynter/vasp.yml.

    Parameters
    ----------
    cfgfile : (str), optional
        Path to configuration file. The default is ~/.pynter/vasp.yml.

    Returns
    -------
    (dict)
        Configuration dictionary.
    """
    if os.path.exists(cfgfile):
        with open(cfgfile,"r") as ymlfile:
            return yaml.load(ymlfile,Loader=yaml.FullLoader) # add Loader to not get warning
    else:
        warnings.warn('%s does not exist. Run "pynter configure" in the terminal to create it. Using default settings for now'%cfgfile)
        return  

# Overwrite initial settings
config = load_config()
if config:
    SETTINGS.update(config)
    
vasp_defaults = load_vasp_default()
if vasp_defaults:
    SETTINGS['vasp'].update(vasp_defaults)



def run_local(cmd,printout=True,dry_run=False,**kwargs):
    """
    Run a command locally with subprocess package.

    Parameters
    ----------
    cmd : (str)
        Command to run.
    printout : (bool), optional
        Print command, output and error. The default is True.
    dry_run : (bool), optional
        Return back the command, without executing it. The default is False.
    **kwargs : (dict)
        Kwargs to pass to subprocess.run().

    Returns
    -------
    stdout
    stderr
    """
    if dry_run:
        if printout:
            print(cmd)
        return cmd, ''
    if printout:
        print("Run command: %s" %cmd)
    command = cmd.split()
    proc = subprocess.run(command, capture_output=True, shell=False, text=True,**kwargs)
    stdout = proc.stdout
    stderr = proc.stderr
    if printout:
        print(stdout)
        if stderr:
            print(stderr)
    return stdout, stderr
