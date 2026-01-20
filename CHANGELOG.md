# metaTAXONx: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 - [2026-01-20]

### `Added`:
    - [#11](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metatonx/-/issues/11) `CONFIGURE` now also checks the classifier compatibility with a dummy sequence file

### `Changed`:
    - [#12](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metatonx/-/issues/12) pipeline now works in the processing given a single `sample_id`.

## v1.1.0 - [2026-01-12]

### `Changed`:
    - [#7](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metatonx/-/issues/7) Moved `git submodule` to it's own github and docker image.
    - [#6](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metatonx/-/issues/6) Moved `biotaviz` to it's own github and docker image.
    - [#8](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metatonx/-/issues/8) Changed repository name from `metataxonomics-DSL2` to `metataxonx`.
    - [#9](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metatonx/-/issues/9) `minmax` uses `biom-format` instead of a custom script.

## v1.0.0 - [2025-11-26]

### `Added`:
    - `CHECK_INPUT` formats reads or samplesheet into meta channel
    - `CONFIGURE` now checks and downloads classifiers
    - paired-end read merging can be done by either `Flash` or `PEAR`
    - Contains `MULTIQC` and `OmicFlow` report generation
    - `Denoise` is it's own subworkflow. Now groovy functions create a mapping file instead of python
    - CI performs a simple test with `-stub-run`

### `Changed`:
    - Removed all python-based modularization and changed it to nextflow

## v1.0.0 - [2025-9-1]

### `Changed`:
    - [#4](https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metatonx/-/issues/4) refactored DSL1 to DSL2 of old metataxonomics pipeline