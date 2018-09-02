# scripts

- `concat_sample_beds.R`
    - concatenate sample `.bed` files
    - remove `chrM` intervals

## binding_dataset

- `bed.py`
- `generate_nonbinding_peaks.py`
    - input
        - BED file with columns chr/start/end/sample_count_distinct/max_MACS_score
        - dm6 chromosome sizes
    - output:
        - BED file with same columns as input, with 2 negative background intervals per original peak.
- `helper_functions.py`
- `peak_features.py`
    - input
        - BED file with columns chr/start/end/sample_count_distinct/max_MACS_score
        - FASTA file (soft-masked) for intervals in BED
        - directory with MAST result BED files (one per motif)
    - output
        - unprocessed data file (see `peaks_merged_features_unprocessed.txt` for details)
- `peak_features_preprocess.py`
    - Subset and split data for training/testing
    - input
        - data file from `macs_features.py`
        - option `--minsamples` (default 1)
            - minimum number of samples for binding intervals to be included (does not apply to non-binding intervals)
        - option `--maxrep` (default 1)
            - maximum proportion of repeats for binding intervals to be included (does not apply to non-binding intervals)
        - subsetting options also remove non-binding intervals that are generated from the binding intervals
        - `--test_size` proportion of data to use for testing set
    - output
        - `--train_out` training data output file
        - `--test_out` testing data output file
- `subset_binding_peaks.py`

## master_table

- `chip_seq.py`
    - `SeqSample` object with sample-specific attributes
- `master_table.py`
    - generate melted master table (`results/ChIP_master_table_samples.txt`) from
        - `etc/MACSscore_summary_valid_merged.bed`
        - `etc/MACSscore_summary_valid_fe.bed`
        - `etc/target.fa`
        - `etc/MACSscore_summary_valid_merged.anno`
- `pivot_master_table.py`
    - generate pivoted master table (`results/ChIP_master_table_{fe/macs}.txt`) from melted format
- `pivot_master_table_source.py`
    - old script to generate master table
- `subset_master_table.py`
    - subset pivoted master table with 3 options:
        - max MACS score column only
        - >1 samples only
        - P53match_score > 6
