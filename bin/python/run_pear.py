#!/usr/bin/env python3

"""
run_pear
---------------
.. module:: run_pear
  :synopsis: Runs pear from command line and generates qiime format mapping file
.. moduleauthor:: Jos Boekhorst 
.. adapted by Alem Gusinac at 02-05-2024

Generate mapping.txt and pear_counts.txt

Typical run::

    run_pear.py -m metadata.tsv

Run the script with '-h' for a list of options.

"""

import glob
import os
import sys
import argparse
from utils.tools import get_file_counts, generate_qiime_mapping, get_metadata, process_write_counts
from multiprocessing import Pool

def write_counts(metadata, outfolder, outname, options):
    # Fetch raw and trimmed files
    input_R1 = glob.glob("*R1*")
    assembled = [os.path.basename(file) for file in glob.glob(f"{outfolder}/*fastq.gz")]
    samples = metadata['SAMPLE-ID'].tolist()

    # Remove failed samples
    samples_cleaned = [sample_name for sample_name in samples for assemble in assembled if assemble.__contains__(sample_name)]
    input_R1_cleaned = [R1 for sample_name in samples_cleaned for R1 in input_R1 if R1.__contains__(sample_name)]

    # sort files
    input_R1_cleaned.sort()
    assembled.sort()
    samples_cleaned.sort()

    # Prepare arguments for parallel processing
    args = []
    for sample_name, input_file, assembled_file in zip(samples_cleaned, input_R1_cleaned, assembled):
        if input_file.startswith(sample_name) and assembled_file.startswith(sample_name):
            args.append((sample_name, input_file, assembled_file, outfolder))

    # Parallelize processing with multiprocessing.Pool
    with Pool(processes=options["cpus"]) as pool:
        results = pool.map(process_write_counts, args)

    # write output
    with open(outname, 'w') as f:
        f.write("\t".join(['sample', 'input', 'assembled']) + '\n')
        for result in results:
            f.write(result)


def run_pear(sample, forward, reverse):
    command = f"pear -f {forward} -r {reverse} -o {sample} -y 25G -j 32 -q 30 -v 35 -p 0.0001 > {sample}_pear_stats.txt"
    result = os.system(command)
    if result != 0:
        sys.stderr.write(f"sample: {sample} failed\n")

def process_pear(args):
    sample, fw, rev = args
    run_pear(sample=sample, forward=fw, reverse=rev)
    os.system(f'mv {sample}.assembled.fastq assembled/{sample}.fastq')
    os.system(f'gzip assembled/{sample}.fastq')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PEAR argument parser", add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m', dest="metadata", help='metadata file', required=True)
    parser.add_argument('--cpu', dest="cpus", help='number of cpu cores', required=True)
    options = vars(parser.parse_args())

    # Creates folder
    os.mkdir("assembled")

    # Fetch metadata, forward and reverse files
    fw_files = glob.glob("*R1*")
    rev_files = glob.glob("*R2*")
    metadata = get_metadata(options['metadata'])
    samples = metadata['SAMPLE-ID'].tolist()

    # sort files accordingly
    fw_files.sort()
    rev_files.sort()
    samples.sort()

    # Run PEAR on sorted samples & files (using sample-id as naming)
    args = [(sample, fw, rev) for (sample, fw, rev) in zip(samples, fw_files, rev_files)]
    
    with Pool(processes=int(options["cpus"])) as pool:
        pool.map(process_pear, args)

        
    # Generate stats and mapping files
    generate_qiime_mapping(metadata=get_metadata(options['metadata']), file_string1="fastq", file_string2="", outdir="assembled", reverse=False)
    write_counts(metadata, 'assembled', 'pear_counts.txt', options)
