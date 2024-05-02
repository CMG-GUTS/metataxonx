#!/usr/bin/env python3
# todoc

"""
run_cutadapt
---------------
.. module:: run_cutadapt
  :synopsis: run cutadapt on paired-end illumina-style data
.. moduleauthor:: Jos Boekhorst
.. Adapted by Alem Gusinac, 02-05-2024

run_cutadapt.py

Typical run::

    run_cutadapt.py -f FWD_PRIMER -r REV_PRIMER -c 64 -m metadata.tsv -s 'paired' --pear 'yes'

Run the script with '-h' for a list of options.

"""

import argparse
import sys, os, glob
from utils.tools import get_metadata, generate_qiime_mapping, get_file_counts

def write_counts(metadata, outfolder, outname, reverse=None):
    # Fetch raw and trimmed files
    untrimmed = glob.glob("*R1*")
    trimmed = [os.path.basename(file) for file in glob.glob(f"{outfolder}/*R1*")]
    samples = metadata['SAMPLE-ID'].tolist()
    if reverse:
        untrimmed += glob.glob("*R2*")
        trimmed += [os.path.basename(file) for file in glob.glob(f"{outfolder}/*R2*")]
        samples += metadata['SAMPLE-ID'].tolist()

    # sort files
    untrimmed.sort()
    trimmed.sort()
    samples.sort()

    # Writes output
    with open(outname, 'w') as f:
        f.write("\t".join(['sample', 'input', 'trimmed']) + '\n')
        for (sample_name, untrim_file, trim_file) in zip(samples, untrimmed, trimmed):
            if untrim_file.startswith(sample_name) and trim_file.startswith(sample_name):
                precount = get_file_counts(untrim_file)
                postcount = get_file_counts(f'{outfolder}/{trim_file}')
                f.write("\t".join([sample_name, str(precount), str(postcount)]) + '\n')

def run_cutadapt(sample, forward, options, reverse=None):
    if reverse:
        command = f"cutadapt -j {options['cpus']} --discard-untrimmed -g {options['fw']} -G {options['rev']} -o trimmed/{sample}_R1.fastq.gz -p trimmed/{sample}_R2.fastq.gz {forward} {reverse}"
    else:
        command = f"cutadapt -j {options['cpus']} --discard-untrimmed -g {options['fw']} -o trimmed/{sample}_R1.fastq.gz {forward}"
    result = os.system(command)
    if result != 0:
        sys.stderr.write("Error running cutadapt, aborting\n")
        sys.exit(1)


############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "run cutadapt on paired-end illumina-style data"


# main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description, add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', dest='fw', help='fwd primer sequence', required=True)
    parser.add_argument('-r', dest='rev', help='rev primer sequence', required=True)
    parser.add_argument('-m', dest="metadata", help='metadata file', required=True)
    parser.add_argument('-s', dest="seq_read", help='single or paired end sequence reading', required=True)
    parser.add_argument('--pear', dest="pear", help="PEAR 'yes' or 'no'", required=True)
    parser.add_argument('-c', dest="cpus", help='number of cpus to us', default=1)
    options = vars(parser.parse_args())

    os.mkdir("trimmed")
    metadata = get_metadata(options['metadata'])
    # Run cutadapt, create qiime mapping file and summarize counts for either paired or single end sequences
    if options['seq_read'] == "paired" or options['pear'] == "yes":
        for (sample, fw, rev) in zip(metadata['SAMPLE-ID'], metadata['FILENAME_fw'], metadata['FILENAME_rev']):
            run_cutadapt(sample=sample, forward=fw, reverse=rev, options=options)
            generate_qiime_mapping(metadata=metadata, file_string1="R1", file_string2="R2", outdir="trimmed", reverse=True)
            write_counts(metadata, "trimmed", "cutadapt_counts.txt", reverse=True)

    elif options['seq_read'] == 'single' and options['pear'] == "no":
        for (sample, fw) in zip(metadata['SAMPLE-ID'], metadata['FILENAME_fw']):
            run_cutadapt(sample=sample, forward=fw, options=options)
            generate_qiime_mapping(metadata=metadata, file_string1="R1", file_string2="R2", outdir="trimmed", reverse=False)
            write_counts(metadata, "trimmed", "cutadapt_counts.txt", reverse=False)
    sys.exit(0)
