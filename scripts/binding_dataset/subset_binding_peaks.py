import argparse
import pandas as pd

from helper_functions import get_genome_seq_records_dict, peak_file_cols


def _subset_binding_peaks(merged_bed, minsamples, rep_cutoff_type,
                          rep_cutoff, output_fp, genome_seq_records):
    df_binding = pd.read_table(merged_bed, header=None,
                               names=('chr', 'start', 'end',
                                      'sample_count_distinct', 'max_MACS_score'))
    # add sequence features
    df_binding['length'] = df_binding['end'] - df_binding['start']
    df_binding['seq_records'] = [genome_seq_records[interval['chr']][interval['start']:interval['end']]
                                 for i, interval in df_binding.iterrows()]
    df_binding['repeat_count'] = [sum(map(seq_record.seq.count, ['g', 'c', 'a', 't']))
                                  for seq_record in df_binding['seq_records']]
    df_binding['repeat_proportion'] = 1. * df_binding['repeat_count'] / df_binding['length']
    # subset by minsamples and repeat proportion cutoff
    if rep_cutoff_type == 'max':
        rows_binding_drop = ((df_binding['sample_count_distinct'] < args.minsamples) |
                             (df_binding['repeat_proportion'] > args.rep_cutoff))
    else:  # min
        rows_binding_drop = ((df_binding['sample_count_distinct'] < args.minsamples) |
                             (df_binding['repeat_proportion'] < args.rep_cutoff))

    df_binding = df_binding[~rows_binding_drop].reset_index(drop=True)
    # drop unnecessary columns
    df_binding.drop(['length', 'seq_records', 'repeat_count', 'repeat_proportion'],
                    axis=1, inplace=True)
    df_binding.insert(3, 'id', df_binding.index)
    assert list(df_binding.columns) == peak_file_cols
    df_binding.to_csv(output_fp, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged_bed', required=True,
                        help='BED file with unlabeled columns chr/start/end/sample_count_distinct/max_MACS')
    parser.add_argument('--minsamples', required=False, type=int, default=1,
                        help='minimum number of samples for binding intervals to be included (does not apply to non-binding intervals)')
    parser.add_argument('--rep_cutoff_type', required=True, type=str, choices=['min', 'max'],
                        help='(min/max) repeat proportion cutoff type for binding intervals to be included')
    parser.add_argument('--rep_cutoff', required=True, type=float,
                        help='cutoff value for rep_cutoff_type')
    parser.add_argument('--genome_fasta', required=True,
                        help='FASTA file for genome (one sequence per chromosome)')
    parser.add_argument('-o', required=True,
                        help='output filepath')
    args = parser.parse_args()

    genome_seq_records = get_genome_seq_records_dict(args.genome_fasta)
    _subset_binding_peaks(args.merged_bed, args.minsamples,
                          args.rep_cutoff_type, args.rep_cutoff, args.o,
                          genome_seq_records)
