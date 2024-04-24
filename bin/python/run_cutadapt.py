#!/usr/bin/env python3
# todoc

"""
run_cutadapt
---------------
.. module:: run_cutadapt
  :synopsis: run cutadapt on paired-end illumina-style data
.. moduleauthor:: Jos Boekhorst

run_cutadapt.py

Typical run::

    run_cutadapt.py -f FWD_PRIMER -r REV_PRIMER -c 64 *.fastq.gz

Run the script with '-h' for a list of options.

"""

import argparse
import sys
import os
import subprocess


def write_counts(infiles, outfolder, outname):
    infiles = [element for element in infiles if "_R2" not in element]
    with open(outname, 'w') as f:
        f.write("\t".join(['sample', 'input', 'trimmed']) + '\n')
        for infile in infiles:
            sample = infile.split("_")[0]
            result = subprocess.run([f'zcat {infile}|wc'], shell=True, capture_output=True, text=True)
            precount = "%i" % (int(result.stdout.split()[0]) / 4)
            result = subprocess.run([f'zcat {outfolder}/{infile}|wc'], shell=True, capture_output=True, text=True)
            postcount = "%i" % (int(result.stdout.split()[0]) / 4)
            f.write("\t".join([sample, str(precount), str(postcount)]) + '\n')


############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "run cutadapt on paired-end illumina-style data"


# main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description, add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', dest='fwd', help='fwd primer sequence', required=True)
    parser.add_argument('-r', dest='rev', help='rev primer sequence', required=True)
    parser.add_argument('-c', dest="cpus", help='number of cpus to us', default=1)
    parser.add_argument('input_files', help="input files", nargs='+')
    options = vars(parser.parse_args())

    os.mkdir("trimmed")
    names = [element.split("_")[0] for element in options['input_files'] if "R1" in element]
    for name in names:
        command = f"cutadapt -j {options['cpus']} --discard-untrimmed -g {options['fwd']} -G {options['rev']} -o trimmed/{name}_R1.fastq.gz -p trimmed/{name}_R2.fastq.gz {name}_R1.fastq.gz {name}_R2.fastq.gz"
        result = os.system(command)
        if result != 0:
            sys.stderr.write("Error running cutadapt, aborting\n")
            sys.exit(1)
    write_counts(options['input_files'], "trimmed", "cutadapt_counts.txt")
    sys.exit()
