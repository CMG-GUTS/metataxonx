#!/usr/bin/env python3
# todoc
"""
validate_mapping
-----------------------

.. module:: validate_mapping
  :synopsis: validates mapping, WUR nextflow edition

Validates mapping file. Checks for:

* presence "sample" column at the start
* duplicate samples
* illegal characters
* empty lines
* presence of fastq files

"""

import string
import sys
import argparse
import os


def load_txt(filename):
    fileinput = open(filename, 'r')
    text = fileinput.read()
    fileinput.close()
    return text

def column_string(n):
    """go from integer to Excel column letter"""
    rest = n
    string = ""
    while rest > 0:
        module = (rest - 1) % 26
        string = chr(65 + module) + string
        rest = int((rest - module) / 26)
    return string


def is_float(errors, cell, column_index, row_index):
    """For "CORRELATION" columns: is it a valid float?"""
    try:
        if cell != "":
            float(cell)
    except ValueError:
        errors.append(
            f'Could not convert to float (not a number): cell {column_string(column_index + 1)}{row_index + 1:d}, "{cell}"\n')
    return errors


def check_cell(errors, cell, allowed_list, column_index, row_index):
    # check for illegal characters in the supplied cell string; used by "check_characters". allowed_list is a list of strings.
    okchars = ""
    for element in allowed_list:
        okchars += element
    for letter in cell:
        if letter not in okchars:
            errors.append(
                'Invalid character in cell %s%i: "%s"\n' % (column_string(column_index + 1), row_index + 1, letter))
    return errors


def check_characters(lines, errors):
    # check for illegal characters
    header = lines[0].split('\t')
    for column_i, column_header in enumerate(header):
        for line_i, line in enumerate(lines):
            if column_header.split('_')[0].upper() in header_keywords:
                cell = line.split('\t')[column_i]
                errors = check_cell(errors, cell, [digits, letters, '_'], column_i, line_i)
    return errors


def check_samples(lines, errors):  # check for duplicate sample names
    samples = [element.split('\t')[0] for element in lines[1:]]
    reported = set([])
    for sample in samples:
        if samples.count(sample) != 1 and sample not in reported:
            errors.append(f'Duplicate sample name: "{sample}"\n')
            reported.add(sample)
    return errors


def check_compulsory_headers(lines, errors):  # check for compulsory first & last columns
    if lines[0].split('\t')[0].upper() not in ['SAMPLE', 'SAMPLE-ID']:
        errors.append(f'Top left cell (A1) must be "Sample" but is ' + lines[0].split("\t")[0] + '\n')
    if lines[0].split('\t')[-1].upper() != 'DESCRIPTION':
        errors.append('Top right cell must be "Description"\n')
    return errors


def get_mapping(column_nr, lines):
    mapping = {}
    for line in lines[1:]:
        lineg = line.split('\t')
        sample = lineg[0]
        data = lineg[column_nr]
        if data != '':
            mapping[sample] = data
    return mapping


def check_files(lines, errors, fastq_path):
    header = lines[0].split('\t')
    for column_i, column_header in enumerate(header):
        for line_i, line in enumerate(lines[1:]):
            if options['dummy']:
            # qiime mode: assume files are put in $PWD by nextflow
                if column_header.split('-')[-1].upper() == "FILEPATH":
                    filename = line.split('\t')[column_i]
                    if not os.path.isfile(filename):
                        errors.append(f"File {filename} does not exist\n")
            else:
                if column_header.split('_')[0].upper() == "FILENAME":
                    filename = fastq_path + '/' + line.split('\t')[column_i]
                    if not os.path.isfile(filename):
                        errors.append(f"File {filename} does not exist\n")
    return errors


def check_duplicate_headers(lines, errors):
    header = lines[0].upper().split('\t')
    original_header = lines[0].split('\t')
    if len(header) != len(set(header)):
        duplicates = []
        for i, element in enumerate(header):
            if header.count(element) != 1:
                duplicates.append(original_header[i])
        duplicates = list(set(duplicates))
        errors.append(f"One or more duplicates in header: {'; '.join(duplicates)}\n")
    return errors


def generate_mapping_qiime(lines):
    with open("mapping_for_nextflow.txt", 'w') as f:
        header = lines[0].split('\t')
        f.write(lines[0] + "\n")
        for line_i, line in enumerate(lines[1:]):
            line = line.split('\t')
            for column_i, column_header in enumerate(header):
                #print(column_header.split('-')[-1].upper())
                cell = line[column_i]
                if column_header.split('-')[-1].upper() == "FILEPATH":
                    cell = "$PWD/" + cell.split('/')[-1]
                    line[column_i] = cell
            f.write('\t'.join(line) + '\n')
    f.close()


################
# main program #
################

# note: default parameters are set in argparse object (top of __main__)
description = "Mapping file validation"
letters = string.ascii_uppercase + string.ascii_lowercase
digits = string.digits
header_keywords = ['RANKSTAT', 'FILE', 'CORRELATION', 'PAIREDGROUPBY', 'PAIREDTIMEPOINT', 'PAIREDCLASS']

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=description, add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', dest='infile', help='name of input file', required=False, default="mapping.txt")
    parser.add_argument('-p', dest='path', help='path to fastq files', required=False, default="fastq")
    parser.add_argument('-d', dest='dummy', help='create dummy filepath columns for qiime', action='store_true', required=False)
    options = vars(parser.parse_args())
    infile = options['infile']
    fastq_path = options['path']

    # read input
    try:
        lines = load_txt(infile).rstrip('\n').split('\n')
    except FileNotFoundError:
        sys.stderr.write(f'ERROR: could not read input file "{infile}"\n')
        sys.exit(1)

    errors = []

    # the actual checks

    errors = check_compulsory_headers(lines, errors)  # check for first & last column (Sample,Description)
    errors = check_characters(lines, errors)  # check for illegal characters
    errors = check_samples(lines, errors)  # check for duplicate samples
    errors = check_files(lines, errors, fastq_path)  # check for exicstence of all input (fastq) files
    errors = check_duplicate_headers(lines, errors)
    # TODO: add validation of the paired rankstat stuff

    # reporting
    if len(errors) == 0:
        sys.stderr.write("Found no errors in mapping file\n")
    else:
        sys.stderr.write("Found %i errors:\n" % len(errors))
        for error in errors:
            sys.stderr.write(error)
        sys.exit(2)

    if options['dummy']:
        # generate a new mapping file for nextflow-qiime
        generate_mapping_qiime(lines)

    sys.exit(0)