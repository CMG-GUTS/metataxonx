[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# Introduction
Running the nextflow pipeline requires the addition of parameters in the nextflow.config. By default the nextflow pipeline deploys singularity images make it possible to run on many clusters of CMBI. 

## Usage
Clone repository:
```bash
git clone --recurse-submodules https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metataxonomics-DSL2.git
```

## Running from command line
```bash
nextflow run metataxonomics-DSL2/qiime.nf -c metataxonomics-DSL2/nextflow.config -work-dir /your/work/dir --cpus 5
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
nextflow run metataxonomics-DSL2/qiime.nf -c metataxonomics-DSL2/nextflow.config -work-dir /your/work/dir --cpus 5 -profile singularity

```