import argparse
import pandas as pd

score_to_colname = {'macs': 'MACS', 'fe': 'FE'}


def pivot_master_table(master_table_filepath, output, score=None):
    """Generate a pivoted master table."""
    peak_cols = ['chr', 'start', 'end',
                 'repeat_count', 'peak_length', 'repeat_proportion',
                 'detailed_annotation', 'annotation', 'gene']
    # TODO: add repeat_type
    mt_df = pd.read_table(master_table_filepath)
    pivot_table = pd.pivot_table(mt_df, index=peak_cols, columns='sample_name',
                                 fill_value=0)
    out_df = pivot_table[score_to_colname[score]].reset_index()
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
    pivot_master_table(args.master_table, args.output, score=args.score)
