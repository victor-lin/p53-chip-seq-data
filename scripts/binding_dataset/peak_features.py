import argparse
import itertools
import numpy as np
import os
import pandas as pd
from collections import OrderedDict
from Bio.SeqUtils import GC

from helper_functions import get_genome_seq_records_dict


reverse_complements = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
k_values = (2, 3, 4, 5, 6)
mast_filename_suffix = '_meme-chip_mast'


def get_df_mast(mast_fp, site_name, df_peaks):
    """Return dataframe with aggregated P53match score and count columns.

    Join with intervals from `df_peaks`, use df_peaks interval identifier as index.

    column : aggregation function
    - P53match_score_max_{site_name} : max
    - P53match_score_sum_{site_name} : sum
    - P53match_count_{site_name} : len
    """
    with open(mast_fp) as f:
        comments, names = itertools.takewhile(lambda line: line.startswith('#'), f)
        names = names[1:].split()

    colname_score_max = 'P53match_score_max_' + site_name
    colname_score_sum = 'P53match_score_sum_' + site_name
    colname_count = 'P53match_count_' + site_name

    df_mast = pd.read_table(mast_fp, comment='#', sep='\s+', index_col=False, header=None, names=names)
    df_mast[colname_score_max] = -np.log10(df_mast['hit_p-value'])
    df_mast[colname_score_sum] = -np.log10(df_mast['hit_p-value'])
    df_mast.rename(columns={'sequence_name': 'chr'}, inplace=True)

    # merge with chr, start, end from df_peaks
    df_merge_mast_all = df_peaks.merge(df_mast, on='chr')
    intersect = (df_merge_mast_all.hit_start >= df_merge_mast_all.start) & (df_merge_mast_all.hit_end <= df_merge_mast_all.end)
    df_merge_mast_intersect = df_merge_mast_all[intersect]
    group_merge_mast_intersect = df_merge_mast_intersect.groupby(list(df_peaks.columns))

    df_mast_cols = group_merge_mast_intersect.agg({'hit_p-value': len,
                                                   colname_score_max: max,
                                                   colname_score_sum: sum}).reset_index()
    df_mast_cols[colname_count] = df_mast_cols['hit_p-value'].astype(int)
    df_mast_cols.set_index(['chr', 'start', 'end'], inplace=True)
    df_mast_cols = df_mast_cols.loc[:, (colname_score_max, colname_score_sum, colname_count)]
    return df_mast_cols


def get_kmer_dict(seq_record, k_values):
    """"Return OrderedDict of k-mer counts.

    Counts are merged forward and backwards strand for reverse complement pairs.
    Remove duplicate kmers after merging reverse complements.

    Parameters
    ----------
    seq_record : Bio.Seq.Seq object
    k_values : iterable
        values of k to determine k-mers
    """
    kmer_dict = OrderedDict()
    for k in k_values:
        kmers = get_kmers_unique(k)
        counts = {kmer: seq_record.seq.upper().count_overlap(kmer) + seq_record.seq.upper().reverse_complement().count_overlap(kmer)
                  for kmer in kmers}
        sum_counts = sum(counts.values())
        if sum_counts == 0:
            kmer_dict.update(OrderedDict((k, 0) for k, v in counts.iteritems()))
        else:
            kmer_dict.update(OrderedDict((k, 1. * v / sum_counts) for k, v in counts.iteritems()))
    return kmer_dict


def get_kmers_unique(n):
    """Return list of unique kmers (no duplicate reverse complements)."""
    kmers_all = [''.join(x) for x in itertools.product('ACGT', repeat=n)]
    kmers_seen = set()
    kmers_unique = []
    for kmer in kmers_all:
        if get_reverse_complement(kmer) in kmers_seen:
            continue
        kmers_seen.add(kmer)
        kmers_unique.append(kmer)
    return kmers_unique


def get_reverse_complement(seq_str):
    """Return reverse complement of a sequence string."""
    return ''.join([reverse_complements[base] for base in seq_str[::-1]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--peaks_all', required=True,
                        help='Text file with columns chr/start/end/id/sample_count_distinct/max_MACS_score/length/repeat_proportion')
    parser.add_argument('--peaks_all_phastcons', required=True,
                        help='BED file with columns chr/start/end/average_phastCon')
    parser.add_argument('--mast_dir', required=True,
                        help='Directory containing MAST files for motifs')
    parser.add_argument('--genome_fasta', required=True,
                        help='FASTA file for genome (one sequence per chromosome)')
    parser.add_argument('-o', required=True,
                        help='Output filepath')
    args = parser.parse_args()

    genome_seq_records = get_genome_seq_records_dict(args.genome_fasta)

    df_all = pd.read_table(args.peaks_all)
    df_all['length'] = df_all['end'] - df_all['start']
    df_all['seq_records'] = [genome_seq_records[interval['chr']][interval['start']:interval['end']]
                             for i, interval in df_all.iterrows()]
    df_all['repeat_count'] = [sum(map(seq_record.seq.count, ['g', 'c', 'a', 't']))
                              for seq_record in df_all['seq_records']]
    df_all['repeat_proportion'] = 1.0 * df_all['repeat_count'] / df_all['length']
    df_all['GC_content'] = [GC(seq_record.seq) for seq_record in df_all['seq_records']]

    print 'adding phastcons...'
    df_phastcons = pd.read_table(args.peaks_all_phastcons, header=None,
                                 na_values='NAN', keep_default_na=False,
                                 names=('chr', 'start', 'end', 'average_phastCon'))
    assert len(df_all) == len(df_phastcons)
    df_all['average_phastCon'] = df_phastcons['average_phastCon']

    print 'adding k-mer features'
    df_kmers = pd.DataFrame([get_kmer_dict(seq_record, k_values)
                             for seq_record in df_all['seq_records']])
    df_all = pd.concat([df_all.reset_index(drop=True), df_kmers], axis=1)

    print 'adding MAST scores...'
    for fn in next(os.walk(args.mast_dir))[2]:
        mast_fp = os.path.join(args.mast_dir, fn)
        site_name = os.path.splitext(fn)[0]
        if mast_filename_suffix in site_name:
            site_name = site_name.split(mast_filename_suffix)[0]
        print '\tadding site %s' % site_name
        df_mast_cols = get_df_mast(mast_fp, site_name, df_all.loc[:, ['chr', 'start', 'end']])
        df_all = df_all.join(df_mast_cols, on=('chr', 'start', 'end'))

    # drop unnecessary columns, fill NA with 0
    df_all.drop(['chr', 'start', 'end', 'seq_records', 'repeat_count'],
                axis=1, inplace=True)
    df_all.fillna(0, inplace=True)
    df_all.to_csv(args.o, sep='\t', index=False)
