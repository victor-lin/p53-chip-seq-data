import argparse
import pandas as pd
import numpy as np

score_to_colname = {'macs': 'MACS', 'fe': 'FE'}


def pivot_master_table(master_table_filepath, output, score):
    """Generate a pivoted master table."""
    peak_cols = ['chr', 'start', 'end',
                 'repeat_count', 'peak_length', 'repeat_proportion',
                 'detailed_annotation', 'annotation', 'repeat_type', 'gene']
    mt_df = pd.read_table(master_table_filepath)
    mt_df['repeat_type'] = mt_df['repeat_type'].fillna('dummy')
    sample_cols = pd.pivot_table(mt_df, index=peak_cols, columns='sample_name',
                                 fill_value=0)
    sample_cols = sample_cols[score_to_colname[score]]
    peak_agg_cols = mt_df.groupby(peak_cols).agg({'MACS': max,
                                                  'FE': max,
                                                  'sample_name': lambda x: x.nunique()})
    peak_agg_cols.rename({'MACS': 'maxMACS',
                          'FE': 'maxFE',
                          'sample_name': 'sample_count'},
                         axis=1, inplace=True)
    out_df = pd.concat([peak_agg_cols, sample_cols], axis=1)
    out_df = out_df.reset_index().replace('dummy', np.nan)
    out_df.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--master_table', required=True,
                        help='Master table filepath')
    parser.add_argument('-o', '--output', required=True,
                        help='Ouptut filepath')
    parser.add_argument('--score', required=True,
                        help='Sample-based score to put under')
    args = parser.parse_args()
    pivot_master_table(args.master_table, args.output, args.score)
