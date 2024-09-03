#!/usr/bin/env python3
# todoc

"""
dummy_metadata
---------------
.. module:: dummy_metadata
  :synopsis: Generate Metadata.tsv from folder with fastq files
.. moduleauthor:: Jos Boekhorst
.. Adapted by Alem Gusinac, 02-05-2024

Generate Metadata.tsv from folder with fastq files

Typical run::

    dummy_metadata.py -s 'paired' --pear 'yes' -o metadata.tsv

Run the script with '-h' for a list of options.

"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys, os, glob, re

############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "Generate Metadata.tsv from folder with fastq files"

# main program
if __name__ == '__main__':
    parser = ArgumentParser(description=description, add_help=True, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', dest='outfile', help='name of output file', default="std")
    parser.add_argument('-f', dest='fullpath', help='use full path for file columns', action='store_true')
    parser.add_argument('-n', dest='nopath', help='leave out path for fiel columns', action='store_true')
    parser.add_argument('-s', dest="seq_read", help='single or paired end sequence reading', required=True)
    parser.add_argument('--pear', dest="pear", help="PEAR 'yes' or 'no'", required=True)
    parser.add_argument('-p', dest='manualpath', help='use supplied path for file columns', default="NA")
    parser.add_argument('-a', dest='single', help='set to no when paired-end data supplied', default="yes")
    parser.add_argument('--fastq', dest='fastq_files', help="fastq files", nargs='+')

    options = vars(parser.parse_args())

    # pattern for fastq files
    #file_extension = r"(\.filt)?\.fastq\.gz"
    file_extension = r"(_R1_)"

    # Fetch forward reads
    fw_files = glob.glob("*R1*")
    # Split only at R1 and remove special characters such as . _ and R1
    samples = [re.split(file_extension, file, maxsplit=1)[0] for file in fw_files]        
    
    # sort forward and samples
    fw_files.sort()
    samples.sort()

    # Write mapping file
    with open(options["outfile"], "w") as outfile:
        # Create header column
        outfile.write("\t".join(['SAMPLE-ID', 'FILENAME_fw', 'FILENAME_rev', 'DESCRIPTION']) + '\n')

        # Paired reads populate both fw and rev columns
        if options['seq_read'] == 'paired' or options['pear'] == 'yes':
            rev_files = glob.glob("*R2*")
            rev_files.sort()

            for (sample, fw, rev) in zip(samples, fw_files, rev_files):
                outfile.write("\t".join([sample, fw, rev, ""]) + '\n')
        
        # Single reads will only populate fw columns
        elif options['seq_read'] == 'single':
            for (sample, fw) in zip(samples, fw_files):
                outfile.write("\t".join([sample, fw, "", ""]) + '\n')
        else:
            print("seq_read should be either 'paired' or 'single'!")
            sys.exit(2)
    sys.exit()