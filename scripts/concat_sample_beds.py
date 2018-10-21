"""
Script to concatenate MACS output files.

ASSUMPTIONS
===========

Each sample has a file in `bed_directory` and `fe_directory`.
- file under `bed_directory` has the filename format `<SAMPLE_NAME>_peaks.bed`
- file under  `fe_directory` has the filename format `<SAMPLE_NAME>_peaks.xls`
- each pair of sample files has the same row order for sequences

Result file has columns
    'chr', 'start', 'end', 'length', 'sample_name', 'MACS', 'FE'
sorted by
    'chr', 'start'
"""

import os
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--bed_directory', required=True,
                    help='BED files directory')
parser.add_argument('-o', '--out_file', required=True,
                    help='output filepath')
parser.add_argument('--fe_directory', required=True,
                    help='FE files directory')
parser.add_argument('--ignore_chr', required=False,
                    help='chromsome to exclude')
args = parser.parse_args()


class NoFEFileError(Exception):
    pass

fnames_bed = os.walk(args.bed_directory).next()[2]
fnames_fe = os.walk(args.fe_directory).next()[2]

assert len(fnames_bed) == len(fnames_fe)

for fname in fnames_bed:
    sample_name = fname.rstrip('_peaks.bed')
    sample_df = pd.read_table(os.path.join(args.bed_directory, fname), names=["chr", "start", "end", "PeakID", "MACS"])
    sample_df.drop('PeakID', axis=1, inplace=True)
    sample_df['sample_name'] = sample_name

    fname_fe = sample_name + '_peaks.xls'
    if fname_fe not in fnames_fe:
        raise NoFEFileError('FEfile {} not found'.format(fname_fe))

    # skip first 23 rows of MACS arguments
    fe_df = pd.read_table(os.path.join(args.fe_directory, fname_fe), skiprows=23)
    assert len(sample_df) == len(fe_df)
    sample_df['FE'] = fe_df['fold_enrichment']

    if fname == fnames_bed[0]:
        summary_df = sample_df.copy()
    else:
        summary_df = summary_df.append(sample_df)

    summary_df['length'] = summary_df['end'] - summary_df['start']
    summary_df = summary_df[['chr', 'start', 'end', 'length', 'sample_name', 'MACS', 'FE']]

if args.ignore_chr:
    summary_df = summary_df.loc[summary_df['chr'] != args.ignore_chr, ]

summary_df.sort_values(by=['chr', 'start'], inplace=True)
summary_df.to_csv(args.out_file, sep='\t', index=False)
