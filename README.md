# p53-chip-seq-data

see wiki for information regarding master table, figures, and files.

Full commands for generating files can be found in [Makefile](Makefile).

## Requirements

### R

* [ggplot2](http://ggplot2.org/)
* [optparse](https://github.com/trevorld/optparse/)

### Python

* [pandas](http://pandas.pydata.org/)
* [numpy](http://www.numpy.org/)
* [pandasql](https://github.com/yhat/pandasql/)
* [Biopython](http://biopython.org/) >= 1.70

### Command line tools

* [bedtools](http://bedtools.readthedocs.io/)
* [BEDOPS](http://bedops.readthedocs.io/)
* [HOMER](http://homer.salk.edu/)

## Binding Dataset

### Quickstart

Generate train/test data files with nonbinding intervals:

    $ make train_test_split

Files are generated in the `etc` directory.

### Data File Format

Files generated when calling `make train_test_split`.

Options (defined in Makefile):

- minimum number of samples (`MINSAMPLES`)
- repeat threshold type (`REP_THRESHOLD_TYPE`)
- repeat threshold cutoff (`REP_CUTOFF`)
- number of nonbinding intervals to generate per binding interval (`N_NONBINDING_INTERVALS`)
	- Searches within [1kb, 10kb] on both sides of the binding interval, then +10kb increments if needed.

Columns:

- binding (0 or 1)
- length
- repeat_proportion
- GC_content
- average_phastCon
- P53_match_count (per motif)
- P53_match_score_max (per motif)
- P53_match_score_sum (per motif)
- 2-mer count proportions (10)
- 3-mer count proportions (32)
- 6-mer count proportions (2080)

### Generated Files

- `etc/peaks_binding_all_samples.txt`
- `etc/peaks_binding_merged_maxMACS.bed`
- `etc/peaks_binding_merged_subset.txt`
- `etc/peaks_nonbinding.txt`
- `etc/peaks_all.txt`
- `etc/peaks_all.bed`
- `etc/peaks_all_phastcons.bed`

- `etc/peaks_merged_features__minsamples_{int}__rep_{max/min}{0-1}__nonbinding_{int}.txt`
- `etc/peaks_merged_features_train__minsamples_{int}__rep_{max/min}{0-1}__nonbinding_{int}.txt`
- `etc/peaks_merged_features_test__minsamples_{int}__rep_{max/min}{0-1}__nonbinding_{int}.txt`

## Master Table

### Quickstart

Generate master table in melted format:

    $ make results/ChIP_master_table_samples.txt

Generate master table in pivoted format:

    $ make results/ChIP_master_table_fe.txt
    $ make results/ChIP_master_table_macs.txt

### Procedure

1. Add FE and peak length columns to sample `.bed` files, concatenate all samples into one file
1. Merge using `bedtools`
    1. `MACSscore_summary_valid_merged.bed`
1. Use merged regions from `MACSscore_summary_valid_merged.bed`
1. Combine information into one table

### Generated Files

- `etc/peaks_binding_all_samples.txt`
- `etc/peaks_binding_merged.bed`
- `etc/peaks_binding_merged.anno`
- `etc/peaks_binding_merged.fa`
