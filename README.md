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

## API

```bash
usage: fastq_handler [-h] -i IN_DIR -o OUT_DIR [-s SLEEP] [-n TAG] [--keep_names] [--monitor] [--max-size MAX_SIZE] [--downsize] [--merge]

Process fastq files.

optional arguments:
  -h, --help            show this help message and exit
  -i IN_DIR, --in_dir IN_DIR
                        Input directory
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory
  -s SLEEP, --sleep SLEEP
                        Sleep time
  -n TAG, --tag TAG     name tag, if given, will be added to the output file names
  --keep_names          keep original file names
  --monitor             run indefinitely
  --max-size MAX_SIZE   max size of the output file, in kilobytes
  --downsize            downsize fastq files to max-size
  --merge               merge files
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

## Acknowledgements

This package took great inspiration inspiration from [fastq-filter](https://github.com/LUMC/fastq-filter), which is a tool for filtering FASTQ files. From it we took the idea of using [xopen](https://github.com/pycompression/xopen) and [dnaio](https://github.com/marcelm/dnaio) in combination as well the actions design pattern. The code for fastq filtering was adapted to the concatenation process.

## Maintainers

- [**@xiaodre21**](https://github.com/xiaodre21)
- [**@santosjgnd**](https://github.com/SantosJGND)
- [**@insaflu**](https://github.com/insapathogenomics)
