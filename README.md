# p53-chip-seq-plots

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

## Input Files

in `data/`

### Sample `.bed` files

`data/ValidSamples/`

* `[sample_name]_peaks.bed`
* Headers absent
* Column order:
 * `Chr`
 * `Start`
 * `End`
 * `PeakID`
 * `MACS_Score`

### FE Files

`data/FEFiles/`

* `[sample_name]_peaks.xls`
* ignore first 23 lines
* column order:
    * `chr`
    * `start` (+1)
    * `stop`
    * `length` (`stop` - `start - 1`)
    * MACS score (col 7)
    * `fold_enrichment`

### `mast_out.bed`
* notable columns:
    * `Seq_name`
    * `hit_start`
    * `hit_end`
    * `hit_p.value`

### `target.fa`

## Procedure

1. Add FE and peak length columns to sample `.bed` files
1. Merge using `bedtools`
    1. `MACSscore_summary_valid_merged.bed`
1. Use merged regions from `MACSscore_summary_valid_merged.bed`
1. `mast_out.bed`: join with merged regions (must be contained)
1. remove ChrM regions

## Generated files

in `etc/`

### `MACSscore_summary_valid_fe.bed`

    concat_sample_beds.R

* concatenation of sample BED files with additional columns for fold enrichment and peak length

### `MACSscore_summary_valid_merged.bed`

    $ bedtools merge -i <(tail -n+2 MACSscore_summary_valid_fe.bed | sort -k1,1 -k2,2n)

* chr, start, stop columns

### `MACSscore_summary_valid_merged.anno`

    $ module load gcc/5.2.0 homer/4.8
    $ annotatePeaks.pl <(cut -f 1,2,3 MACSscore_summary_valid_merged.bed) dm6


## Master Tables

generated in `results/`

### `master_table_macs.txt`

* **chr**/**start**/**end** from merged regions (`MACSscore_summary_valid_merged.bed`)
* Column for each sample name, with respective **MACS Score** under the column if the sample has an intersection with the merged region
* **Detailed Annotation** / **Gene Name** from annotation of merged regions (`MACSscore_summary_valid_merged.anno`)
* **Max Matrix Score** from intersections with `mast_out.bed`
* **Max MACS Score** from intersections with samples
* **Max Fold Enrichment** from intersections with  samples
* **Repeats**: proportion of repeats from `target.fa`

### `master_table_fe.txt`
* same as `master_table_macs.txt` but with **Fold Enrichment** under sample name columns

## Figures

Figures can be generated using `data_plots.R`:

    > source("data_plots.R")

There will be a file prompt. Choose any master_table file.

Use any pre-defined plot:

1. `plot_chr_vs_macs(rep_cut)`
1. `plot_fe_vs_matrix()`
1. `plot_macs_vs_fe()`
1. `plot_macs_vs_matrix()`
1. `plot_matrix_score_hist(rep_cut)`
1. `plot_rep_cutoffs(matrix_cut, interval)`

Or create your own using the base plot:

    plot_base() + [ggplot layers]
