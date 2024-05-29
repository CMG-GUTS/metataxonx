import pandas as pd
import os, sys, glob, subprocess

def get_file_counts(filepath):
    """
    Calls Unix command 'wc' and counts the number of lines

    Parameters:
        - filepath (str)    : absolute path to file

    Returns:
        - count (Int)       : Integer number    

    """
    result = subprocess.run([f'zcat {filepath} | wc'], shell=True, capture_output=True, text=True)
    count = "%i" % (int(result.stdout.split()[0]) / 4) 
    return count

def get_metadata(infile):
    """
    Loads metadata into a pandas dataframe

    Parameters:
        - infile (str)      : Filename

    Returns:
        - metadata (Pandas.DataFrame) : Tabular structure    

    """
    if os.path.split(infile)[1].split('.')[-1] == "csv":
        metadata = pd.read_csv(infile, sep=",")
    elif os.path.split(infile)[1].split('.')[-1] == "tsv":
        metadata = pd.read_csv(infile, sep="\t")
    else:
        raise IOError("File format is not accepted; metadata should be .tsv or .csv file format")
    return metadata

def add_path(df, path_name, read_end):
    """
    Rewrites pandas dataframe column by their absolute filepath

    Parameters:
        - df (Pandas.DataFrame) : Standard Metataxonomics table in the format of at least: 
            'SAMPLE-ID' 'FILENAME_fw' 'FILENAME_rev' 'DESCRIPTION'
        - path_name (str)       : Path, such as current directory to location of the filename

    Returns:
        - df (Pandas.DataFrame) : modified dataframe  

    """
    for i in df.index:
        file_fw = path_name + "/" + df.loc[i, 'FILENAME_fw']
        df.loc[i, 'FILENAME_fw'] = file_fw

        if read_end == "paired":
            file_rev = path_name + "/" + df.loc[i, 'FILENAME_rev']
            df.loc[i, 'FILENAME_rev'] = file_rev
    return df

def update_metadata(metadata, sample_id, file_path, column_name, outdir):
    """
    Checks the existence of a file, and updates it's new filepath

    Parameters:
        - metadata (Pandas.DataFrame) : Standard Metataxonomics table in the format of at least: 
            'SAMPLE-ID' 'FILENAME_fw' 'FILENAME_rev' 'DESCRIPTION'
        - sample_id (str)   : sample-id of the table, default = 'SAMPLE-ID'
        - file_path (str)   : path to the filename
        - column_name (str) : Name of column to be modified
        - outdir (str)      : Folder containing the output files

    Returns:
        - metadata (Pandas.DataFrame) : modified dataframe  

    """
    for i in metadata.index:
        if metadata.loc[i, sample_id] in file_path:
            metadata.loc[i, column_name] = os.getcwd() + f"/{outdir}/" + file_path
    return metadata

def generate_qiime_mapping(metadata, outdir, file_string1, file_string2, reverse=None):
    """
    Generates a qiime mapping file for either single or paired-end reads

    Parameters:
        - metadata (Pandas.DataFrame)   : Standard Metataxonomics table in the format of at least: 
            'SAMPLE-ID' 'FILENAME_fw' 'FILENAME_rev' 'DESCRIPTION'
        - outdir (str)                  : Folder containing the output files
        - file_string1 (str)            : String to recognize forward filename
        - file_string2 (str)            : String to recognize reverse filename
        - reverse (bool)                : If reverse is not None or not False, paired-end qiime mapping is generated, otherwise single-end

    Returns:
        - Creates tab-separated 'mapping.txt' file for qiime  

    """
    for files in os.listdir(f"{outdir}/"):
        if file_string1 in files:
            metadata = update_metadata(metadata=metadata, sample_id='SAMPLE-ID', file_path=files, column_name='FILENAME_fw', outdir=outdir)
        if reverse and file_string2 in files:
            metadata = update_metadata(metadata=metadata, sample_id='SAMPLE-ID', file_path=files, column_name='FILENAME_rev', outdir=outdir)

    columns = ['SAMPLE-ID', 'FILENAME_fw', 'DESCRIPTION']
    if reverse:
        columns.append('FILENAME_rev')
        df_subset = metadata[columns].rename(columns={'FILENAME_fw' : 'forward-absolute-filepath', 'FILENAME_rev' : 'reverse-absolute-filepath'})
        df_final = df_subset[['SAMPLE-ID', 'forward-absolute-filepath', 'reverse-absolute-filepath', 'DESCRIPTION']]
    else:
        df_subset = metadata[columns].rename(columns={'FILENAME_fw' : 'absolute-filepath'})
        df_final = df_subset[['SAMPLE-ID', 'absolute-filepath', 'DESCRIPTION']]
    
    # saving metadata files
    df_final.to_csv("mapping.txt", sep='\t', index=False, header=True)

def checksZeroDivision(num1, num2):
    if num1 == 0.0 or num2 == 0.0:
        return 0.0
    else:
        return num1 / num2