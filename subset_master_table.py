#!/usr/bin/env python
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--master_table', required=True,
                    help='Master table filepath (file in pivoted format)')
parser.add_argument('--max_out', required=False,
                    help='Master table subset: max scores only')
parser.add_argument('--intersect_out', required=False,
                    help='Master table subset: peaks with >1 samples only')
parser.add_argument('--six_mat_out', required=False,
                    help='Master table subset: peaks with max P53match_score>6')
args = parser.parse_args()

mt_df = pd.read_table(args.master_table)
non_max_cols = [col for col in mt_df.columns if '_all' not in col]
if args.max_out:
    mt_df[non_max_cols].to_csv(args.max_out, sep='\t', index=False)
if args.intersect_out:
    mt_df.loc[mt_df['P53match_score_max'] > 6,
              non_max_cols].to_csv(args.intersect_out, sep='\t', index=False)
if args.six_mat_out:
    mt_df.loc[mt_df['P53match_score_max'] > 6,
              non_max_cols].to_csv(args.intersect_out, sep='\t', index=False)
