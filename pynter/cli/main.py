#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:59:39 2023

@author: villa
"""

import argparse
import os
import warnings

from pynter.cli.automations import setup_automation
from pynter.cli.analysis import setup_analysis
from pynter.cli.config import setup_config, run_config
from pynter.cli.defects import setup_defects
from pynter.cli.hpc import setup_hpc
from pynter.cli.inputs import setup_inputs
from pynter.cli.materials_project import setup_mp_database
from pynter.cli.plotter import setup_plotter
from pynter.cli.slurm import setup_job_script
from pynter.cli.phase_diagram import setup_phase_diagram


def main():
    """
    Handle main.
    """    
    parser = argparse.ArgumentParser(
        description="""
        This script works based on several sub-commands with their own options.
        To see the options for the sub-commands, type "pynter sub-command -h".
        """
        )

    subparsers = parser.add_subparsers()
    
    homedir = os.getenv("HOME")
    config_exists = os.path.exists(os.path.join(homedir,'.pynter/config.yml'))
    config_vasp_exists = os.path.exists(os.path.join(homedir,'.pynter/vasp.yml'))
    if config_exists is False and config_vasp_exists is False:
        warnings.warn('configuration files do not exist, run "pynter configure" in the terminal to create them')


    setup_analysis(subparsers)
    setup_automation(subparsers)  
    setup_config(subparsers)
    setup_defects(subparsers)
    setup_hpc(subparsers)
    setup_inputs(subparsers)
    setup_job_script(subparsers)
    setup_mp_database(subparsers)
    setup_phase_diagram(subparsers)
    setup_plotter(subparsers)
    
    args = parser.parse_args()
    
    try:
        args.func
    except AttributeError:
        parser.print_help()
        raise SystemExit("Please specify a command.")
    return args.func(args)

if __name__ == '__main__':
    main()