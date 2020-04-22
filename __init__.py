import yaml


def load_config(cfgfile="config.yml"):
    with open(cfgfile,"r") as ymlfile:
        return yaml.load(ymlfile)