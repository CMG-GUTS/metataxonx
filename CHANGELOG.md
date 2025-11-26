# metaTAXONx: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2025-11-26]

### `Added`
    - `CONFIGURE` now checks and downloads classifiers
    - paired-end read merging can be done by either `Flash` or `PEAR`
    - Contains `MULTIQC` and `OmicFlow` report generation
    - `Denoise` is it's own subworkflow. Now groovy functions create a mapping file instead of python

### `Changed`:
    - Removed all python-based modularization and changed it to nextflow
