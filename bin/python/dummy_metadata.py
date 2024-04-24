#!/usr/bin/env python3
# todoc

"""
dummy_metadata
---------------
.. module:: dummy_metadata
  :synopsis: Generate Metadata.tsv from folder with fastq files
.. moduleauthor:: Jos Boekhorst

Generate Metadata.tsv from folder with fastq files

Typical run::

    dummy_metadata.py -o Metadata.tsv fastq/*

Run the script with '-h' for a list of options.

"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import os

############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "Generate Metadata.tsv from folder with fastq files"

illegal_characters = [".", "_"]

# main program
if __name__ == '__main__':
    parser = ArgumentParser(description=description, add_help=True, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', dest='outfile', help='name of output file', default="std")
    parser.add_argument('-f', dest='fullpath', help='use full path for file columns', action='store_true')
    parser.add_argument('-n', dest='nopath', help='leave out path for fiel columns', action='store_true')
    parser.add_argument('-p', dest='manualpath', help='use supplied path for file columns', default="NA")
    parser.add_argument('-a', dest='single', help='set to no when paired-end data supplied', default="yes")
    parser.add_argument('fastq_files', help="fastq files", nargs='+')

    options = vars(parser.parse_args())

    names = []
    infiles = options['fastq_files']
    infiles.sort()
    samples = set([])

    postfixF = ""
    if infiles[0].find("_R1.") != -1:
        postfixF = infiles[0][infiles[0].find("_R1."):]
    if infiles[0].find(".R1.") != -1:
        postfixF = infiles[0][infiles[0].find(".R1."):]
    if infiles[0].find("_1.") != -1:
        postfixF = infiles[0][infiles[0].find("_1."):]
    if infiles[0].find("R1_") != -1:
        sys.stderr.write("Detected BaseClear-style run labels (sample_R1_001.xxx), rename with JOS_rename_baseclear.py and try again. Aborting\n")
        sys.exit(1)
    if options['single'] == "no" and postfixF == "":
        sys.stderr.write("Could not determine postfix, looked for R1. and _1.\nAborting\n")
        sys.exit(1)
    if options['single'] == "yes" and postfixF == "":
        postfixF = "NA" # nanopore data: singletons
    if options['single'] == "no":
        text = ['\t'.join(['Sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath', 'DESCRIPTION'])]
    else:
        text = ['\t'.join(['Sample-id', 'absolute-filepath', 'DESCRIPTION'])]
    for name in infiles:
        if options['fullpath']:
            current = os.getcwd() + '/'
            name = current + name.replace(current, "")
        nice_name = name.split('/')[-1].split('_')[0]
        if name.find(postfixF) != -1 or postfixF == "NA":
            for character in illegal_characters:
                if postfixF == "NA":
                    nice_name = nice_name.split('.')[0]
                if nice_name.replace(postfixF, "").count(character) != 0:
                    sys.stderr.write(f"Invalid character {character} in {name}, aborting\n")
                    sys.exit(1)
            fwd_file = name
            if options['single'] == "no":
                rev_file = fwd_file.replace(postfixF, postfixF.replace('1', '2'))
                if not os.path.isfile(rev_file):
                    sys.stderr.write(f"Could not find reverse file going with {name}, aborting\n")
                    sys.exit(1)
                if options['manualpath'] != "NA":
                    fwd_file = options['manualpath'] + "/" + fwd_file.split('/')[-1]
                    rev_file = options['manualpath'] + "/" + rev_file.split('/')[-1]
                nice_name = name.split('/')[-1].split('_')[0]
                text.append('\t'.join([nice_name, fwd_file, rev_file, nice_name]))
            elif options['single'] == "yes":
                if options['manualpath'] != "NA":
                    fwd_file = options['manualpath'] + "/" + fwd_file.split('/')[-1]
                nice_name = name.split('/')[-1].split('_')[0]
                filename = fwd_file.replace(postfixF.split('.')[0], "")
                text.append('\t'.join([nice_name, filename, nice_name]))

    text = '\n'.join(text)
    if options['outfile'] == 'std':
        print(text)
    else:
        with open(options['outfile'], 'w') as f:
            f.write(text + '\n')
    sys.exit()
