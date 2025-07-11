[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://code.askimed.com/nf-test)

## Introduction **metaTAXONx**

The 16S pipeline is a best-practice suite for the preprocessing, denoising and annotation of sequencing data obtained via 16S rRNA marker-gene sequencing. The pipeline contains [NF-core modules](https://github.com/nf-core/modules) and other local modules that are in the similar format. It can be runned via both docker and singularity containers.

## Pipeline summary
The pipeline is able to perform different taxonomic annotation on either (single/paired) reads. The different subworkflows can be defined via `--bypass_<method>` flags, a full overview is shown by running `--help`.

The pipeline performs preprocessing of the reads via the removal of primers or adapters via [cutadapt](https://cutadapt.readthedocs.io/en/stable/) and paired-end read merging via [PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html).Before and after each step the quality control will be assessed via [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and a [multiqc](https://github.com/MultiQC/MultiQC) report is created as output. The denoising of single-end reads is performed via [DADA2](https://benjjneb.github.io/dada2/) in batches or in paralell with the module [run-dada2-batch](https://github.com/agusinac/run-dada2-batch).

**Taxonomy assignment**

**Diversity analysis**

**Report/Visualisation**


## Installation
> [!NOTE]
> Make sure you have installed the latest [nextflow](https://www.nextflow.io/docs/latest/install.html#install-nextflow) version! 

Clone the repository in a directory of your choice:
```bash
git clone https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metataxonomics-DSL2.git
```

The pipeline is containerised, meaning it can be runned via docker or singularity images. No further actions need to be performed when using the docker profile, except a docker registery needs to be set on your local system, see [docker](https://docs.docker.com/engine/install/). In case singularity is used, please specify the `singularity.cacheDir` in the `nextflow.config` so that singularity images are saved there and re-used again.

## Usage
Since the latest version, metaBIOMx works with both a samplesheet (CSV) format or a path to the input files. Preferably, samplesheets should be provided.
```bash
nextflow run main.nf --input <samplesheet.csv> -work-dir work -profile singularity
nextflow run main.nf --input <'*_{1,R1,2,R2}.{fq,fq.gz,fastq,fastq.gz}'> -work-dir work -profile singularity
```

## Support

If you are having issues, please [create an issue](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metataxonomics-DSL2/-/issues)

## Citations

You can cite the `metataxonx` using the following DOI: 

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md)
file.