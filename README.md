# Introduction
Running the nextflow pipeline requires the addition of parameters in the nextflow.config. By default the nextflow pipeline deploys singularity images make it possible to run on many clusters of CMBI. 

## Running from command line
```bash
nextflow run nextflow/metataxonomics-DSL2/qiime.nf -c nextflow/metataxonomics-DSL2/nextflow.config -work-dir /your/work/dir --cpus 5
```

## Running as a jobscript
Make sure to have the cores specified in ``-pe smp`` the same as in ``--cpus``.
In case nextflow is in your anaconda/miniconda environment, it is required to use ``eval "$(conda shell.bash hook)"`` prior to ``conda activate`` and ``nextflow run``.
```bash
#!/bin/bash
#$ -N jobname
#$ -o jobname.out
#$ -e jobname.err
#$ -V
#$ -q all.q@narrativum.umcn.nl
#$ -pe smp 5

eval "$(conda shell.bash hook)" 
conda activate nextflow
cd $HOME
nextflow run nextflow/metataxonomics-DSL2/qiime.nf -c nextflow/metataxonomics-DSL2/nextflow.config -work-dir /your/work/dir --cpus 5

```