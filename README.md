# p53-chip-seq-data

see wiki for information regarding master table, figures, and files.

Full commands for generating files can be found in [Makefile](Makefile).

## Requirements

### R

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

Generate train/test data files with non-binding intervals:

    $ make train_test_split

### Data File Format

Files generated when calling `make train_test_split`.

Options (defined in Makefile):

- minimum number of samples (`MINSAMPLES`)
- repeat threshold type (`REP_THRESHOLD_TYPE`)
- repeat threshold cutoff (`REP_CUTOFF`)
- number of non-binding intervals to generate on each side per binding interval (`N_NONBINDING_INTERVALS`)
	- Searches within [1kb, 10kb] on both sides of the binding interval, then +10kb increments if needed.

Columns:

- binding (0 or 1)
- length
- repeat_proportion
- GC_content
- average_phastCon
- P53match_count (per motif)
- P53match_score_max (per motif)
- P53match_score_sum (per motif)
- 2-mer count proportions (10)
- 3-mer count proportions (32)
- 6-mer count proportions (2080)

### Generated Files

1. Full unprocessed data file: `etc/peaks_merged_features__minsamples_{int}__rep_{max/min}{0-1}__nonbinding_{int}.txt`
2. Training set: `results/datafiles/peaks_merged_features__minsamples_{int}__rep_{max/min}{0-1}__nonbinding_{int}-preprocessed_train.txt`
3. Testing set: `results/datafiles/peaks_merged_features__minsamples_{int}__rep_{max/min}{0-1}__nonbinding_{int}-preprocessed_test.txt`

The full file (1) is unprocessed.

The `TEST_SIZE` variable in `Makefile` determines train/test split ratio for files (2) and (3). Default 1:1. Split is based on binding intervals only – non-binding intervals follow the same split as the respective binding interval it was derived from. These files are used in the machine learning process described below.

### Machine Learning Analysis

Model used: [`sklearn.svm.SVC`](http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html) – Support Vector Classifier, based on Support Vector Machine.

Ranking mthod: [`sklearn.feature_selection.RFE`](http://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.RFE.html) – Recursive Feature Elimination.

## Master Table

### Quickstart

Generate master table with either max FE or max MACS score under sample columns:

    $ make results/datafiles/ChIP_peaks_master_table_fe.txt
    $ make results/datafiles/ChIP_peaks_master_table_macs.txt

### Procedure

1. Add FE and peak length columns to sample `.bed` files, concatenate all samples into one file
1. Merge using `bedtools`
    1. `MACSscore_summary_valid_merged.bed`
1. Use merged regions from `MACSscore_summary_valid_merged.bed`
1. Combine information into one table
