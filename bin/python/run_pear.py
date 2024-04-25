#!/usr/bin/env python3
"""minimal Python wrapper for Pear to prevent doing stuff in nextflow Java"""
import glob
import os
import sys
import subprocess


def write_counts(infiles, outfolder, outname):
    with open(outname, 'w') as f:
        f.write("\t".join(['sample', 'input', 'assembled']) + '\n')
        for infile in infiles:
            sample = infile.split("_")[0]
            result = subprocess.run([f'zcat {infile}|wc'], shell=True, capture_output=True, text=True)
            precount = "%i" % (int(result.stdout.split()[0]) / 4)
            result = subprocess.run([f'zcat {outfolder}/{infile.replace("_R1", "")}|wc'], shell=True, capture_output=True, text=True)
            postcount = "%i" % (int(result.stdout.split()[0]) / 4)
            f.write("\t".join([sample, str(precount), str(postcount)]) + '\n')


if __name__ == "__main__":
    all_files = glob.glob("*fastq*")
    os.mkdir("assembled")
    for name in all_files:
        fwd = ""
        rev = ""
        if name.find("_R1.") != -1:
            fwd = name
            rev = name.replace("_R1.", "_R2.")
        elif name.find("_1.") != -1:
            fwd = name
            rev = name.replace("_1.", "_2.")
        nicename = '_'.join(name.split('.')[0].split('_')[:-1])
        if fwd != "":
            command = f"""pear \
            -f {fwd} \
            -r {rev} \
            -o {nicename} \
            -y 25G \
            -j 32 \
            -q 30 \
            -v 35 \
            -p 0.0001 \
            > {nicename}_pear_stats.txt"""
            print(command)
            result = os.system(command)
            if result != 0:
                sys.stderr.write("Error running pear, aborting\n")
                sys.exit(2)
            os.system(f'mv {nicename}.assembled.fastq assembled/{nicename}.fastq')
            os.system(f'gzip assembled/{nicename}.fastq')
    write_counts(glob.glob("*R1.fastq.gz"), 'assembled', 'pear_counts.txt')
