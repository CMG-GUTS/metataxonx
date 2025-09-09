[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
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
git clone --recursive https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metataxonomics-DSL2.git
```

In case you want to update a current git clone with all it's submodules please use
```bash
git submodule update --recursive --remote
```

The pipeline is containerised, meaning it can be runned via docker or singularity images. No further actions need to be performed when using the docker profile, except a docker registery needs to be set on your local system, see [docker](https://docs.docker.com/engine/install/). In case singularity is used, images are automatically cached within the project directory.

## Usage
Since the latest version, metaBIOMx works with both a samplesheet (CSV) format or a path to the input files. Preferably, samplesheets should be provided.
```bash
nextflow run main.nf --input <samplesheet.csv> -work-dir work -profile singularity
nextflow run main.nf --input <'*_{1,R1,2,R2}.{fq,fq.gz,fastq,fastq.gz}'> -work-dir work -profile singularity
```

### 📋 Sample Metadata File Specification

metaTAXONx expects your sample input data to follow a **simple, but strict** structure to ensure compatibility and allow upfront validation. The input should be provided as a **CSV** file where **each entry = one sample** with specified sequencing file paths. Additional properties not mentioned here will be ignored by the validation step.

---

### **Minimum requirement**
- **`sample_id`** ➡ every entry **must** have a unique, non-empty sample identifier.
- No spaces are allowed in sample IDs — use underscores `_` or dashes `-` instead.
- **`forward_read`** ➡ every entry **must** provide a path to an existing forward read FASTQ file (gzipped).
- If `reverse_read` is provided, `forward_read` must also be present.
Example:

| sample_id | forward_read | reverse_read |
|-----------|---------------|--------------------|
| sample1   | sample1_R1.fastq.gz | sample1_R2.fastq.gz |
| sample_2  | D029327_1.fastq.gz | D029327_2.fastq.gz |
| S3        | L9283_R1.fastq.gz | L9283_R1.fastq.gz |

---

### **Properties and Validation Rules**

#### 🔹 Required properties

| Property     | Type   | Rules / Description                                                                                   |
|--------------|--------|----------------------------------------------------------------------------------------------------|
| `sample_id`     | string | Unique sample ID with no spaces (`^\S+$`). Serves as an identifier.                                  |
| `forward_read` | string | File path to forward sequencing read. Must be non-empty string matching FASTQ gzipped pattern. File must exist. |

#### 🔹 Optional property

| Property       | Type   | Rules / Description                                                                                   |
|----------------|--------|----------------------------------------------------------------------------------------------------|
| `reverse_read` | string | File path to reverse sequencing read. Same constraints as `forward_read`. Required if specified.   |

## Support

If you are having issues, please [create an issue](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metataxonomics-DSL2/-/issues)

## Citations

You can cite the `metataxonx` using the following DOI: 

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md)
file.