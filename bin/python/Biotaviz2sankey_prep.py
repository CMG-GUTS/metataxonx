#!/usr/bin/env python3
# todoc
"""
Biotaviz2sankey_python
---------------
.. module:: Biotaviz2sankey_prep
  :synopsis: Make input for R sankey diagram from biotaviz-style compositional table
.. moduleauthor:: Jos Boekhorst

Make sankey diagram from biotaviz-style compositional table

Typical run::

    scriptname -i BiotaViz_relative_abundance.txt

Run the script with '-h' for a list of options.

"""
import sys
import argparse

description = "Make input for R sankey diagram from biotaviz-style compositional table"


def read_biotaviz(filename):
    """read biotaviz-style compositional table"""
    traces = []
    data = {}
    with open(filename, "r") as f:
        lines = f.read().rstrip().split('\n')
    samples = lines[0].split('\t')[2:]
    for line in lines[1:]:
        abundances = {}
        lineg = line.split('\t')
        trace = lineg[0]
        if trace[-1] != '.':
            trace += '.'  # Excel removes the final . in a cell when it could be a number (as in 1.)
        label = lineg[1]
        scores = lineg[2:]
        for i, score in enumerate(scores):
            abundances[samples[i]] = float(scores[i])
        data[trace] = [label, abundances]
        traces.append(trace)
    return data, traces, samples


def filter_traces(data, traces, dataset, threshold):
    """remove traces of taxa not meeting relative abundance threshold"""
    filtered_traces = []
    average_abundances = {}
    for trace in traces:
        abundances = []
        for sample in dataset:
            abundances.append(data[trace][1][sample])
        average_abundance = sum(abundances) / float(len(abundances))
        if average_abundance >= threshold:
            average_abundances[trace] = average_abundance
            filtered_traces.append(trace)
    return filtered_traces, average_abundances


def print_csv(data, traces, outfile_prefix, datasets, threshold):
    """generate comma-separated sankey network file"""
    for dataset in datasets:
        filtered_traces, average_abundances = filter_traces(data, traces, datasets[dataset], threshold)
        outfile = f"{outfile_prefix}_{dataset}.csv"
        with open(outfile, "w") as f:
            f.write(",".join(["link1", "link2", "label", "value"]) + '\n')
            for i, trace in enumerate(filtered_traces):
                label, abundances = data[trace]
                abundance = average_abundances[trace]
                label = label.replace(' - ', '~').split("~")[1]
                parent = '.'.join(trace.split('.')[:-2])+'.'
                try:
                    parent_index = filtered_traces.index(parent)
                except ValueError:  # it was the root node, so no parent
                    parent_index = -1
                    link1 = "link1"
                    link2 = "link2"
                if parent_index != -1:
                    link1 = parent_index
                    link2 = i
                f.write(",".join([f"{link1}", f"{link2}", f" {label}:{abundance*100:.1f}%", f"{abundance:.4f}"]) + '\n')


def get_sets(filename):
    """parse mapping file: what subsets should be compared"""
    with open(filename, "r") as f:
        lines = f.read().rstrip().split('\n')
    headers = {}
    datasets = {}
    for i, element in enumerate(lines[0].split('\t')):
        if element.split('_')[0].lower() == 'rankstat':
            dataset = '_'.join(element.split("_")[1:])
            headers[i] = dataset
            datasets[dataset] = {}
    for line in lines[1:]:
        lineg = line.split('\t')
        sample_id = lineg[0]
        for i, element in enumerate(lineg):
            if i in list(headers.keys()) and element != '':
                sample_class = element
                dataset = headers[i]
                if sample_class not in datasets[dataset]:
                    datasets[dataset][sample_class] = []
                datasets[dataset][sample_class].append(sample_id)
    flattened = {}
    for setname in datasets:  # as sankey is per-class, the second level complicates scripting and is removed
        for classname in datasets[setname]:
            flattened[setname + '_' + classname] = datasets[setname][classname]
    return flattened


# main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-i', dest='infile', help='name of input file', required=True)
    parser.add_argument('-o', dest='outfile_prefix', help='prefix for output files', required=True)
    parser.add_argument('-m', dest='mappingfile', help='name of mapping file', required=False, default="")
    parser.add_argument('-t', dest='min_abundance', help='abundance threshold', required=False, default=0.01, type=float)

    options = vars(parser.parse_args())
    data, traces, samples = read_biotaviz(options['infile'])

    if options['mappingfile'] != "":
        datasets = get_sets(options['mappingfile'])
    else:
        datasets = {}
    datasets['All'] = samples  # alwas have a set containing all samples (removes requirement of always having a mapping file with an all column)
    print_csv(data, traces, options['outfile_prefix'], datasets, options['min_abundance'])
    
    sys.exit(0)
    
