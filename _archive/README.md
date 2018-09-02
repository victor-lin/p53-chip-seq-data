# _archive

- `concat_mast.sh`
    - concatenate all files in `data/MAST`, write to `etc/mast_concat.bed`. Add column for sample name (derived from filename)
- `concat_sample_annos.py`
    - concatenate sample `.anno` files
    - remove `chrM` intervals
- `data_plots.R`
    - pre-defined plotting functions, given a master table file
