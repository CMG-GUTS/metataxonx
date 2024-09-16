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
import sys, os, glob, re
from utils.tools import get_metadata, generate_qiime_mapping, get_file_counts
from multiprocessing import Pool

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

    print("Samples:", samples)
    print("Untrimmed files:", untrimmed)
    print("Trimmed files:", trimmed)

    # Writes output
    with open(outname, 'w') as f:
        f.write("\t".join(['sample', 'input', 'trimmed']) + '\n')
        for (sample_name, untrim_file, trim_file) in zip(samples, untrimmed, trimmed):
            if untrim_file.startswith(sample_name) and trim_file.startswith(sample_name):
                precount = get_file_counts(untrim_file)
                postcount = get_file_counts(f'{outfolder}/{trim_file}')
                f.write("\t".join([sample_name, str(precount), str(postcount)]) + '\n')


def adapter_trimming(sample, forward, options, reverse=None):
    if os.path.basename(options['fw_adapter']) == "adapters.fasta":
        adapter_option = f"--discard-untrimmed --minimum-length 20 -a file:{options['fw_adapter']}"
        if reverse:
            adapter_option += f" -A file:{options['rev_adapter']}"
    elif re.search("[ATGC]", options['fw_adapter']):
        adapter_option = f"--discard-untrimmed --minimum-length 20 -a {options['fw_adapter']}"
        if reverse:
            adapter_option += f" -A {options['rev_adapter']}"
    else:
        adapter_option = ""

    command = f"cutadapt -j 4 {adapter_option} -o trimmed/{sample}_R1.fastq.gz -p trimmed/{sample}_R2.fastq.gz {forward} {reverse}"
    return(command)  

def primer_trimming(sample, forward, options, reverse=None):
    if os.path.basename(options['fw_primer']) == "primers.fasta":
        primer_option = f"--discard-untrimmed --minimum-length 20 -g file:{options['fw_primer']}"
        if reverse:
            primer_option += f" -G file:{options['rev_primer']}"
    elif re.search("[ATGC]", options['fw_primer']):
        primer_option = f"--discard-untrimmed --minimum-length 20 -g {options['fw_primer']}"
        if reverse:
            primer_option += f" -G {options['rev_primer']}"
    else:
        primer_option = ""

    command = f"cutadapt -j 4 {primer_option} -o trimmed/{sample}_R1.fastq.gz -p trimmed/{sample}_R2.fastq.gz {forward} {reverse}"
    return(command)

def run_cutadapt(sample, forward, options, reverse=None):
    if options['fw_adapter'] != "skip" and options['fw_primer'] != "skip":
        adapter_cmd = adapter_trimming(sample, forward, options, reverse)
        primer_cmd = primer_trimming(sample, "./trimmed/" + sample + "_R1.fastq.gz", options, "./trimmed/" + sample + "_R2.fastq.gz")
        result = os.system(adapter_cmd)
        if result != 0:
            sys.stderr.write("Error running cutadapt, aborting\n")
            sys.exit(1)
        result = os.system(primer_cmd)   
    elif options['fw_adapter'] != "skip":
        cmd = adapter_trimming(sample, forward, options, reverse)
        result = os.system(cmd)
    elif options['fw_primer'] != "skip":
        cmd = primer_trimming(sample, forward, options, reverse)
        result = os.system(cmd)
    else:
        # place holder returns untrimmed
        cmd = primer_trimming(sample, forward, options, reverse)
        result = os.system(cmd)
    if result != 0:
        sys.stderr.write("Error running cutadapt, aborting\n")
        sys.exit(1)

def process_sample(args):
    sample, fw, rev, options, metadata = args
    run_cutadapt(sample=sample, forward=fw, reverse=rev, options=options)
    generate_qiime_mapping(metadata=metadata, file_string1="R1", file_string2="R2", outdir="trimmed", reverse=rev)

############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "run cutadapt on paired-end illumina-style data"


# main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description, add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-fp', dest='fw_primer', help='fwd primer sequence', required=True)
    parser.add_argument('-rp', dest='rev_primer', help='rev primer sequence', required=True)
    parser.add_argument('-fa', dest='fw_adapter', help='fwd adapter sequence', required=True)
    parser.add_argument('-ra', dest='rev_adapter', help='rev adapter sequence', required=True)
    parser.add_argument('-m', dest="metadata", help='metadata file', required=True)
    parser.add_argument('-s', dest="seq_read", help='single or paired end sequence reading', required=True)
    parser.add_argument('--pear', dest="pear", help="PEAR 'yes' or 'no'", required=True)
    parser.add_argument('-c', dest="cpus", help='number of cpus to us', default=1)
    options = vars(parser.parse_args())

    os.mkdir("trimmed")
    metadata = get_metadata(options['metadata'])
    if options['seq_read'] == "paired" or options['pear'] == "yes":
        samples = zip(metadata['SAMPLE-ID'], metadata['FILENAME_fw'], metadata['FILENAME_rev'])
        args = []
        for sample, fw, rev in samples:
            args.append((sample, fw, rev, options, metadata))

        with Pool(processes=int(options["cpus"])) as pool:
            pool.map(process_sample, args)
        
        # Create stats output
        write_counts(metadata, "trimmed", "cutadapt_counts.txt", reverse=True)
    
    elif options['seq_read'] == 'single' and options['pear'] == "no":
        samples = zip(metadata['SAMPLE-ID'], metadata['FILENAME_fw'])
        args = []
        for sample, fw in samples:
            args.append((sample, fw, None, options, metadata))

        with Pool(processes=int(options["cpus"])) as pool:
            pool.map(process_sample, args)

        # Create stats output
        write_counts(metadata, "trimmed", "cutadapt_counts.txt", reverse=False)


