#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:59:39 2023

@author: villa
"""

import argparse

from pynter.cli.automations import setup_automation
from pynter.cli.analysis import setup_analyse_vasprun
from pynter.cli.plotter import setup_plot_dos, setup_plot_dos_bs, setup_plot_bs, setup_plot_neb
from pynter.cli.slurm import setup_job_script


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
    
    setup_analyse_vasprun(subparsers)
    
    setup_plot_bs(subparsers)
    setup_plot_dos(subparsers)
    setup_plot_dos_bs(subparsers)
    setup_plot_neb(subparsers)
    
    setup_job_script(subparsers)
    
    args = parser.parse_args()
    
    try:
        args.func
    except AttributeError:
        parser.print_help()
        raise SystemExit("Please specify a command.")
    return args.func(args)

if __name__ == '__main__':
    main()