#!/usr/bin/env python3
# todoc

"""
merge_readstats
---------------
.. module:: merge_readstats
  :synopsis: merge parallel_dada2, pear, cutadapt read stats
.. moduleauthor:: Jos Boekhorst
.. adapted by Alem Gusinac at 03-10-2024

Add Pear-assembler and cutadapt statistics to the dada2_denoise stats file of parallel_dada2

Typical run::

    merge_readstats -d denoising_stats.txt -s paired

Run the script with '-h' for a list of options.
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd

############
# SETTINGS #
############
description = "Add Pear-assembler and cutadapt statistics to the dada2_denoise stats file of qiime2"
single_headers = ["sample", "input", "trimmed", "trimmed%", "filtered", "filtered%", "denoised", "denoised%", "nonchim", "nonchim%"]
paired_headers = ["sample", "input", "trimmed", "trimmed%", "assembled", "assembled%", "filtered", "filtered%", "denoised", "denoised%", "nonchim", "nonchim%"]

# main program
if __name__ == '__main__':
    parser = ArgumentParser(description=description, add_help=True,
                                     formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', dest='denoisefile', help='denoise_stats.txt file', default=False)
    parser.add_argument('-c', dest='cutadaptfile', help='cutadapt counts file', default=False)
    parser.add_argument('-p', dest='pearfile', help='pear counts file', default=False)
    parser.add_argument('-o', dest='outfile', help='output filename', default=False)
    parser.add_argument('-s', dest="seq_read", help="paired or single -end reads", required=True)

    options = vars(parser.parse_args())

    # Fetch files if present
    if options["cutadaptfile"]:
        cutadapt_stats = pd.read_csv(options["cutadaptfile"], sep = "\t")

    if options["denoisefile"]:
        denoise_stats = pd.read_csv(options["denoisefile"], sep = "\t")
        denoise_stats["sample"] = denoise_stats.index
        denoise_stats = denoise_stats.drop(columns = ["input"])

    if options["pearfile"]:
        pear_stats = pd.read_csv(options["pearfile"], sep = "\t")
        pear_stats = pear_stats.drop(columns = ["input"])

    # Combine stats for single-end seq
    if options["seq_read"] == "single":
        # merge single stats
        single_stats = cutadapt_stats.merge(denoise_stats, on = "sample")
        baseline = single_stats["input"]

        # Compute percentage differences
        single_stats["trimmed%"] = single_stats["trimmed"] / baseline * 100
        single_stats["filtered%"] = single_stats["filtered"] / baseline * 100
        single_stats["denoised%"] = single_stats["denoised"] / baseline * 100
        single_stats["nonchim%"] = single_stats["nonchim"] / baseline * 100

        # rearranges column order
        single_stats_final = single_stats[single_headers]

        # outputs tsv file
        single_stats_final.to_csv("read_stats.tsv", sep = "\t", index = False)

    if options["seq_read"] == "paired":
        paired_stats = cutadapt_stats.merge(pear_stats, on = "sample").merge(denoise_stats, on = "sample")
        baseline = paired_stats["input"]

        # Compute percentage differences
        paired_stats["trimmed%"] = paired_stats["trimmed"] / baseline * 100
        paired_stats["assembled%"] = paired_stats["assembled"] / baseline * 100
        paired_stats["filtered%"] = paired_stats["filtered"] / baseline * 100
        paired_stats["denoised%"] = paired_stats["denoised"] / baseline * 100
        paired_stats["nonchim%"] = paired_stats["nonchim"] / baseline * 100

        # rearranges column order
        paired_stats_final = paired_stats[paired_headers]

        # outputs tsv file
        paired_stats_final.to_csv("read_stats.tsv", sep = "\t", index = False)