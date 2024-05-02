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
import pandas as pd
from pandas.api.types import is_numeric_dtype, is_string_dtype, is_float_dtype
from utils.tools import get_metadata, add_path

def check_filename_reads(df, errors, seq_arg, pear_arg):
    for col in df.keys():
        if col.lower().find("FILENAME_fw") == -1 and not df["FILENAME_fw"].isnull().any():
            pass
        else:
            errors.append("metadata filename_fw column missing or empty!")
        if seq_arg == "paired" or pear_arg == "yes":
            if col.lower().find("FILENAME_rev") == -1 and not df["FILENAME_rev"].isnull().any():
                pass 
            else:
                errors.append("metadata filename_rev column missing or empty!")
    return errors

def check_compulsory_headers(df, errors, seq_read, pear):
    # Checks compulsory columns; SAMPLE and Description
    if df.keys()[0].lower() != "sample-id" and df.keys()[-1].lower() != "description":
        errors.append("First and Last column should be sample-id and description!")
    
    # Check if single or paired reads are in correct format
    errors = check_filename_reads(df, errors, seq_arg=seq_read, pear_arg=pear)
    
    return df, errors

def check_columns(df):
    # Check existence of header keywords and valid values, corrects if it's not the case
    # TODO: Implement letters & digits check for each column
    for col in df.keys():
        if col.split("_")[0].upper() == "RANKSTAT":
            if not is_string_dtype(df[col]):
                df[col] = df[col].astype(str)
        elif col.split("_")[0].upper() == "CORRELATION":
            if not is_float_dtype(df[col]) and not is_numeric_dtype(df[col]):
                df[col] = df[col].astype(float)
        elif col.split("_")[0].upper() == "PAIREDGROUPBY":
            if not is_string_dtype(df[col]):
                df[col] = df[col].astype(str)
        elif col.split("_")[0].upper() == "PAIREDTIMEPOINT":
            if not is_string_dtype(df[col]):
                df[col] = df[col].astype(str)
        elif col.split("_")[0].upper() == "PAIREDCLASS":
            if not is_string_dtype(df[col]):
                df[col] = df[col].astype(str)
    return df

def generate_corrected_metadata(df, errors, seq_read):
    # Collects cleaned metadata
    df_clean = check_columns(df) 
    # Drops duplicates
    df_clean.drop_duplicates()

    # Create clean metadata format
    df_clean_metadata = add_path(df=df_clean.copy(), path_name=os.getcwd(), read_end=seq_read)
    df_clean_metadata.to_csv("metadata_clean.tsv", sep='\t', index=False, header=True)

    return errors

################
# main program #
################

# note: default parameters are set in argparse object (top of __main__)
description = "Mapping file validation"
letters = string.ascii_uppercase + string.ascii_lowercase
digits = string.digits

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description, add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', dest='infile', help='name of input file', required=True, default="mapping.txt")
    parser.add_argument('-p', dest='path', help='path to fastq files', required=False, default="fastq")
    parser.add_argument('-s', dest="seq_read", help='single or paired end sequence reading', required=True)
    parser.add_argument('--pear', dest="pear", help="PEAR 'yes' or 'no'", required=True)
    parser.add_argument('-d', dest='dummy', help='create dummy filepath columns for qiime', action='store_true', required=False)
    options = vars(parser.parse_args())
    infile = options['infile']

    # read input
    try:
        metadata = get_metadata(infile)
    except FileNotFoundError:
        sys.stderr.write(f'ERROR: could not read input file "{infile}"\n')
        sys.exit(1)

    errors = []

    # the actual checks
    checked_metadata, errors = check_compulsory_headers(metadata, errors, options['seq_read'], options['pear']) # checks for sample-id, filename_fw, filename_rev
    errors = generate_corrected_metadata(checked_metadata, errors, options['seq_read'])

    # reporting
    if len(errors) == 0:
        sys.stderr.write("Found no errors in mapping file\n")
    else:
        sys.stderr.write("Found %i errors:\n" % len(errors))
        for error in errors:
            sys.stderr.write(error)
        sys.exit(2)
    sys.exit(0)