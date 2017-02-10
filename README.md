# p53-chip-seq-plots

see [wiki](https://github.com/zhoulab/p53-chip-seq-plots/wiki) for information regarding master table, figures, and files.

## Requirements

### R

* [ggplot2](http://ggplot2.org/)
* [optparse](https://github.com/trevorld/optparse/)

### Python

* [pandas](http://pandas.pydata.org/)
* [numpy](http://www.numpy.org/)
* [pandasql](https://github.com/yhat/pandasql/)
* [Biopython](http://biopython.org/)

### Command Line

* [bedtools](http://bedtools.readthedocs.io/)
* [Homer](http://homer.salk.edu/)

## Quickstart

    $ module load homer bedtools gcc python R
    $ source /ufrc/zhou/share/virtualenvs/ve/bin/activate

Generate master table in melted format:

    $ make results/ChIP_master_table.txt

Ignore samples: 

    $ make results/ChIP_master_table.txt ignore=<sample1>,<sample2>

Generate master table in pivoted format:

    $ make results/ChIP_master_table_fe.txt
    $ make results/ChIP_master_table_macs.txt

## Procedure

1. Add FE and peak length columns to sample `.bed` files
1. Merge using `bedtools`
    1. `MACSscore_summary_valid_merged.bed`
1. Use merged regions from `MACSscore_summary_valid_merged.bed`
1. `mast_out.bed`: join with merged regions (must be contained)
1. remove ChrM regions
