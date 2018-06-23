import argparse
import sys
import pandas as pd

from helper_functions import get_genome_seq_records_dict, peak_file_cols


def _subset_binding_peaks(merged_bed, minsamples, maxrep, genome_fasta):
    genome_seq_records = get_genome_seq_records_dict(genome_fasta)

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
    # drop unnecessary columns
    # subset by minsamples and maxrep
    rows_binding_drop = ((df_binding['sample_count_distinct'] < args.minsamples) |
                         (df_binding['repeat_proportion'] > args.maxrep))

    df_binding = df_binding[~rows_binding_drop].reset_index(drop=True)
    df_binding.drop(['length', 'seq_records', 'repeat_count', 'repeat_proportion'],
                    axis=1, inplace=True)
    df_binding.insert(3, 'id', df_binding.index)
    assert list(df_binding.columns) == peak_file_cols
    df_binding.to_csv(sys.stdout, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged_bed', required=True,
                        help='BED file with unlabeled columns chr/start/end/sample_count_distinct/max_MACS')
    parser.add_argument('--minsamples', required=False, type=int, default=1,
                        help='minimum number of samples for binding intervals to be included (does not apply to non-binding intervals)')
    parser.add_argument('--maxrep', required=False, type=float, default=1.0,
                        help='maximum proportion of repeats for binding intervals to be included')
    parser.add_argument('--genome_fasta', required=True,
                        help='FASTA file for genome (one sequence per chromosome)')
    args = parser.parse_args()

    _subset_binding_peaks(args.merged_bed, args.minsamples, args.maxrep, args.genome_fasta)
