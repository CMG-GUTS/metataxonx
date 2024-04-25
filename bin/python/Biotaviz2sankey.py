#!/usr/bin/python3
# todoc

"""
Biotaviz2sankey
---------------
.. module:: Biotaviz2sankey
  :synopsis: Create sundquist-style sankey diagram from BiotaViz file
.. moduleauthor:: wrapper Jos, actual scripts by Harm

Typical run::

    scriptname -i Biotaviz_relative.txt -m mapping.txt

Run the script with '-h' for a list of options.

"""

import argparse
import sys
import os


def run_command(command):
    sys.stderr.write(f"Executing: {command}\n")
    os.system(command)


############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "Create sundquist-style sankey diagram from BiotaViz file"

# main program
if __name__ == '__main__':
    work_path = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-i', dest='infile', help='name of input file', required=True)
    parser.add_argument('-m', dest='mapping', help='name of mapping file', required=True)
    options = vars(parser.parse_args())

    run_command(f"python3 {work_path}/sankey-file-prep.py 0.01 false {options['mapping']} false {options['infile']}")
    run_command(
        f"Rscript {work_path}/../R/sankey-diagram-html-generator.R biotaviz_sankey_prepfile-AverageAllSamples.csv")
    run_command(
        f"Rscript {work_path}/../R/sankey-diagram-png-generator.R biotaviz_sankey_prepfile-AverageAllSamples.html")
    sys.exit()
