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

    $ make results/ChIP_master_table_samples.txt

Generate master table in pivoted format:

    $ make results/ChIP_master_table_fe.txt
    $ make results/ChIP_master_table_macs.txt

## Procedure

1. Add FE and peak length columns to sample `.bed` files, concatenate all samples into one file
1. Merge using `bedtools`
    1. `MACSscore_summary_valid_merged.bed`
1. Use merged regions from `MACSscore_summary_valid_merged.bed`
1. Combine information into one table

## scripts

`concat_sample_beds.R`

- concatenate sample `.bed` files
- remove `chrM` intervals

`concat_sample_annos.py`

- concatenate sample `.anno` files
- remove `chrM` intervals

`master_table.py`

- generate melted master table (`results/ChIP_master_table_samples.txt`) from
    - `etc/MACSscore_summary_valid_merged.bed`
    - `etc/MACSscore_summary_valid_fe.bed`
    - `etc/target.fa`
    - `etc/MACSscore_summary_valid_merged.anno`

`pivot_master_table.py`

- generate pivoted master table (`results/ChIP_master_table_{fe/macs}.txt`) from melted format

`mast_concat.sh`

- concatenate all files in `data/MAST`, write to `etc/mast_concat.bed`. Add column for sample name (derived from filename)
