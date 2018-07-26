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

## Quickstart

Generate train/test data files with nonbinding intervals:

    $ make train_test_split

Files are generated in the `etc` directory.

## Data File Format

Files generated when calling `make train_test_split`.

Options (defined in Makefile):

- minimum number of samples
- repeat threshold type
- repeat threshold cutoff
- number of nonbinding intervals to generate per binding interval (within 10kb)

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

## Generated Files

These files are written to the `etc` directory. Columns for non-header BED files are described in parentheses.

- `peaks_all_samples.bed`
    - un-merged ChIP-seq intervals from all samples (`chr/start/end/length/sample_name/MACS_score/FE_score`)
- `peaks_merged_3col.bed`
    - merged ChIP-seq intervals from all samples (`chr/start/end`)
- `peaks_merged_3col_with_negative.bed`
    - merged ChIP-seq intervals from all samples and randomly generated non-binding intervals (`chr/start/end`)
- `peaks_merged_maxMACS.bed`
    - merged ChIP-seq intervals from all samples with MACS score aggregated by max value (`chr/start/end/sample_count_distinct/max_MACS_score`)
- `peaks_merged_maxMACS_with_negative.bed`
    - merged ChIP-seq intervals from all samples and randomly generated non-binding intervals. Non-binding intervals have `max_MACS_score = 0`. (`chr/start/end/sample_count_distinct/max_MACS_score`)
- `peaks_merged.anno`
    - annotation file from HOMER suite's `annotatePeaks.pl`
- `peaks_merged_mask.fa`
    - repeat-masked FASTA corresponding to intervals in `peaks_merged_3col.bed`
- `peaks_merged_softmask_negative.fa`
    - soft-masked FASTA corresponding to intervals in `peaks_merged_3col.bed`
- `peaks_phastcons_with_negative.bed`
    - `peaks_merged_3col_with_negative.bed` with additional column for average phastCon score. Intervals with missing score values have `NAN` values (`chr/start/end/average_phastCon_score`)
- `peaks_merged_features_unprocessed.txt`
    - unprocessed feature file with columns for MACS score (class label), length, repeat proportion, GC content, k-mers (2,3,6), average conservation score, P53 match count and score for each motif