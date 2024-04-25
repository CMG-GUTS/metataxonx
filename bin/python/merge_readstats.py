#!/usr/bin/env python3
# todoc

"""
merge_readstats
---------------
.. module:: merge_readstats
  :synopsis: merge Qiime2 nextflow read accounting data
.. moduleauthor:: Jos Boekhorst

Add Pear-assembler and cutadapt statistics to the dada2_denoise stats file of qiime2

Typical run::

    merge_readstats -d denoising_stats.txt -i fastq/

Run the script with '-h' for a list of options.
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def read_stats(filename, relevants):
    """read counts from txt file"""
    with open(filename, 'r') as f:
        lines = f.read().rstrip().split('\n')
    lines = [element for element in lines if element[0] != "#"]  # skip comments
    header = lines[0].split('\t')
    data = {}
    for line in lines[1:]:  # header skipped
        lineg = line.split('\t')
        sample = lineg[0]
        data[sample] = {}
        for element in relevants:
            try:
                i = header.index(element)
                data[sample][element] = lineg[i]
            except ValueError:  # header not in this file
                pass
    return data


def add_fractions(data, max_header):
    """after each header add percentage of input"""
    new_data = {}
    final_header = []
    for sample in data:
        final_header = []
        new_data[sample] = {}
        new_data[sample]["input"] = data[sample]['input']
        for column in max_header[1:]:  # skip the "input"
            if column in data[sample]:
                final_header.append(column)
                new_data[sample][column] = data[sample][column]
                final_header.append(f"{column}%")
                new_data[sample][f"{column}%"] = "%.2f" % (
                            float(data[sample][column]) / float(data[sample]['input']) * 100)
    return new_data, final_header


def write_stats(filename, new_data, final_header):
    with open(filename, 'w') as f:
        f.write("\t".join(["sample"] + final_header) + '\n')
        for sample in new_data:
            f.write("\t".join([sample] + [new_data[sample][element] for element in final_header]) + '\n')


def combine(denoise_data, cutadapt_data=None, pear_data=None):
    combined_data = denoise_data
    if cutadapt_data is not None:
        for sample in cutadapt_data:
            combined_data[sample]['input'] = cutadapt_data[sample]['input']
            combined_data[sample]['trimmed'] = cutadapt_data[sample]['trimmed']
    if pear_data is not None:
        if cutadapt_data is None:
            for sample in pear_data:
                combined_data[sample]['input'] = pear_data[sample]['input']
        for sample in pear_data:
            combined_data[sample]['assembled'] = pear_data[sample]['assembled']
    return combined_data


############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "Add Pear-assembler and cutadapt statistics to the dada2_denoise stats file of qiime2"
starting_header_relevants = ["sample", "input", "denoised", "merged", "non-chimeric", "trimmed", "assembled"]
max_header = ["sample", "input", "trimmed", "assembled", "denoised", "merged", "non-chimeric"]

# main program
if __name__ == '__main__':
    parser = ArgumentParser(description=description, add_help=True,
                                     formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', dest='denoisefile', help='denoise_stats.txt file', default="denoising_stats.txt")
    parser.add_argument('-c', dest='cutadaptfile', help='cutadapt counts file', default="NA")
    parser.add_argument('-p', dest='pearfile', help='pear counts file', default="NA")
    parser.add_argument('-o', dest='outfile', help='output filename', default="counts.txt")

    options = vars(parser.parse_args())

    denoise_data = read_stats(options['denoisefile'], starting_header_relevants)
    if options['cutadaptfile'] != "NA":
        cutadapt_data = read_stats(options['cutadaptfile'], starting_header_relevants)
    else:
        cutadapt_data = None
    if options['pearfile'] != "NA":
        pear_data = read_stats(options['pearfile'], starting_header_relevants)
    else:
        pear_data = None
    combined_data = combine(denoise_data, cutadapt_data=cutadapt_data, pear_data=pear_data)
    final_data, final_header = add_fractions(combined_data, max_header)
    write_stats(options['outfile'], final_data, final_header)
