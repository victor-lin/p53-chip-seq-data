# p53-chip-seq-data

see wiki for information regarding master table, figures, and files.

## Requirements

### R

* [ggplot2](http://ggplot2.org/)
* [optparse](https://github.com/trevorld/optparse/)

### Python

* [pandas](http://pandas.pydata.org/)
* [numpy](http://www.numpy.org/)
* [pandasql](https://github.com/yhat/pandasql/)
* [Biopython](http://biopython.org/) >= 1.70

### Command Line

* [bedtools](http://bedtools.readthedocs.io/)
* [Homer](http://homer.salk.edu/)

## Quickstart

Generate MACS data file with negatives:

    $ make results/ChIP_MACS_merged_features.txt

Generate master table in pivoted format:

    $ make results/ChIP_master_table_fe.txt
    $ make results/ChIP_master_table_macs.txt

## Result Files

### `ChIP_MACS_merged_features.txt`

Data file with columns:

- max_MACS_score
- length
- repeat_proportion
- GC_content
- 2-mer count proportions (10)
- 3-mer count proportions (32)
- P53_match_count (per motif)
- P53_match_score (per motif)

For each MACS-detected peak from ChIP-seq data, two negative background intervals are generated within 10kb of the original peak.

## Scripts

### `concat_sample_beds.R`

- concatenate sample `.bed` files
- remove `chrM` intervals

### `concat_sample_annos.py`

- concatenate sample `.anno` files
- remove `chrM` intervals

### `master_table.py`

- generate melted master table (`results/ChIP_master_table_samples.txt`) from
    - `etc/MACSscore_summary_valid_merged.bed`
    - `etc/MACSscore_summary_valid_fe.bed`
    - `etc/target.fa`
    - `etc/MACSscore_summary_valid_merged.anno`

### `pivot_master_table.py`

- generate pivoted master table (`results/ChIP_master_table_{fe/macs}.txt`) from melted format

### `concat_mast.sh`

- concatenate all files in `data/MAST`, write to `etc/mast_concat.bed`. Add column for sample name (derived from filename)

### `macs_features.py`

- input
    - BED file with columns chr/start/end/sample_count_distinct/max_MACS_score
    - FASTA file (soft-masked) for intervals in BED
    - directory with MAST result BED files (one per motif)
- output
    - data file (see `ChIP_MACS_merged_features.txt` for details)

### `generate_negative_intervals.py`

- input
    - BED file with columns chr/start/end/sample_count_distinct/max_MACS_score
    - dm6 chromosome sizes
- output:
    - BED file with same columns as input, with 2 negative background intervals per original peak.