import argparse
import csv
import sys
import numpy as np
import pandas as pd

rs = 0  # random seed


def output_df_with_negatives(merged_bed_macs, chrom_sizes):
    df_peaks = pd.read_table(merged_bed_macs, header=None,
                             names=['chr', 'start', 'end', 'sample_count_distinct', 'max_MACS_score'])
    df_peaks['length'] = df_peaks['end'] - df_peaks['start']

    # chromosome size constraints define low and high range for random interval shift.
    with open(chrom_sizes) as f:
        chrom_sizes_dict = {row[0]: int(row[1])
                            for row in csv.reader(f, delimiter='\t')}
    np.random.seed(rs)

    # random interval from left
    df_negative1 = df_peaks.copy()
    df_negative1['sample_count_distinct'] = 0
    df_negative1['max_MACS_score'] = 0
    df_negative1['high'] = df_negative1.apply(lambda row: min(row['start'] - row['length'],
                                                              10000), axis=1)
    df_negative1['low'] = 1000
    df_negative1['rand_shift'] = df_negative1.apply(lambda row: np.random.randint(1000, row['high']) if row['low'] < row['high'] else np.nan, axis=1)
    df_negative1['start'] = df_peaks['start'] - df_negative1['rand_shift'] - df_peaks['length']
    df_negative1['end'] = df_peaks['start'] - df_negative1['rand_shift']
    df_negative1.drop(['high', 'low', 'rand_shift'], axis=1, inplace=True)
    df_negative1.dropna(inplace=True)

    # random interval from right
    df_negative2 = df_peaks.copy()
    df_negative2['sample_count_distinct'] = 0
    df_negative2['max_MACS_score'] = 0
    df_negative2['high'] = df_negative2.apply(lambda row: min(chrom_sizes_dict[row['chr']] - row['end'] - row['length'],
                                                              10000), axis=1)
    df_negative2['low'] = 1000
    df_negative2['rand_shift'] = df_negative2.apply(lambda row: np.random.randint(row['low'], row['high']) if row['low'] < row['high'] else np.nan, axis=1)
    df_negative2['start'] = df_peaks['end'] + df_negative2['rand_shift']
    df_negative2['end'] = df_peaks['end'] + df_negative2['rand_shift'] + df_peaks['length']
    df_negative2.drop(['high', 'low', 'rand_shift'], axis=1, inplace=True)
    df_negative2.dropna(inplace=True)

    df_peaks = df_peaks.append(df_negative1)
    df_peaks = df_peaks.append(df_negative2)

    df_peaks.drop('length', axis=1, inplace=True)
    df_peaks.to_csv(sys.stdout, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--macs_bed', required=True,
                        help='BED file with columns chr/start/end/sample_count_distinct/max_MACS')
    parser.add_argument('--chrom_sizes', required=True,
                        help='Soft-masked FASTA file for intervals in merged_bed')
    args = parser.parse_args()
    output_df_with_negatives(args.macs_bed, args.chrom_sizes)
