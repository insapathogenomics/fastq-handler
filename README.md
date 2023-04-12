# Fastq Handler

[![PyPI version](https://badge.fury.io/py/fastq-handler.svg)](https://badge.fury.io/py/fastq-handler)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/fastq-handler.svg)](https://pypi.python.org/pypi/fastq-handler/)
[![PyPI license](https://img.shields.io/pypi/l/fastq-handler.svg)](https://pypi.python.org/pypi/fastq-handler/)
[![PyPI format](https://img.shields.io/pypi/format/fastq-handler.svg)](https://pypi.python.org/pypi/fastq-handler/)

_A python module to process ONT fastq files by concatenating reads as they are generated during a sequencing run_

## INTRODUCTION

The _fastqc-handler_ module screens folders and subfolders for fastq (fastq or fastq.gz format) files and concatenates them iteratively. This is useful for merging same-sample reads that are split into multiple files, as commonly obtained in ONT sequencing. The output is a one file per
output fastq.gz, containing all reads from the previous files. The output directory structure is maintained.

## INPUT

A directory containing fastq files. The files can be in subfolders (each representing a different sample). The files can be gzipped or not.

## USAGE

```bash
usage: fastq_handler [-h] [-i INPUT] [-o OUTPUT] [-n TAG] [--keep_names]

parse arguments

optional arguments:
    -h, --help            show this help message and exit
    -i INPUT, --input INPUT
                        Input directory
    -o OUTPUT, --output OUTPUT
                        Output directory
    -n TAG, --tag TAG     Tag to add to output file name
    --keep_names          Keep original file names in output file
    --max-size MAX_SIZE   max size of the output file, in kilobytes
```

## REQUIREMENTS

**Modules**

- dataclasses==0.6
- natsort==8.3.1
- pandas==1.5.3
- setuptools==57.4.0
- xopen==1.7.0

## INSTALLATION

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install fastq-handler
```

## MAIN OUTPUTS

> **Note:** The output directory structure is maintained.

- **fastq.gz** files containing all reads from the previous files.
- **log.txt** file containing the concatenation process.

## Maintainers

- [**@xiaodre21**](https://github.com/xiaodre21)
- [**@santosjgnd**](https://github.com/SantosJGND)
- [**@insaflu**](https://github.com/insapathogenomics)
