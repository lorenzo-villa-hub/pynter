#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:59:39 2023

@author: villa
"""

import argparse

from pynter.cli.automations import setup_automation
from pynter.cli.analysis import setup_analysis
from pynter.cli.defects import setup_defects
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
    
    setup_automation(subparsers)
    setup_analysis(subparsers)
    setup_defects(subparsers)
    setup_inputs(subparsers)
    setup_mp_database(subparsers)
    setup_plotter(subparsers)
    setup_job_script(subparsers)
    setup_phase_diagram(subparsers)
    
    args = parser.parse_args()
    
    try:
        args.func
    except AttributeError:
        parser.print_help()
        raise SystemExit("Please specify a command.")
    return args.func(args)

if __name__ == '__main__':
    main()