
import yaml
import os


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
        raise FileNotFoundError('%s does not exist. Run "pynter configure" in the terminal to create it.'%cfgfile)
        return        