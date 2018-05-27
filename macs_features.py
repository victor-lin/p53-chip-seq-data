import argparse
import os
import itertools
import pandas as pd
import numpy as np
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqUtils import GC


reverse_complements = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}


def output_features(merged_bed, merged_fasta_softmask, mast_dir, output_fp):
    """Output data file

    Columns
    - max_MACS_score
    - length
    - repeat_proportion
    - GC_content
    - conservation score
    - kmer_counts


    Parameters
    ----------
    merged_bed : str
        filepath of BED file with columns chr/start/end/sample_count_distinct/max_MACS_score
    merged_fasta_mask : str
        filepath of masked FASTA file for intervals in merged_bed
    merged_fasta_softmask : str
        filepath of soft-masked FASTA file for intervals in merged_bed
    mast_dir : str
        path of directory containing MAST files for motifs
    output_fp : str
        filepath of output file
    """
    # generic features
    df_peaks = pd.read_table(merged_bed, header=None,
                             names=('chr', 'start', 'end', 'sample_count_distinct', 'max_MACS_score'))
    df_peaks['length'] = df_peaks['end'] - df_peaks['start']
    df_peaks['seq_records'] = list(SeqIO.parse(merged_fasta_softmask, 'fasta'))

    df_peaks['repeat_count'] = [sum(map(seq_record.seq.count, ['g', 'c', 'a', 't']))
                                for seq_record in df_peaks['seq_records']]
    df_peaks['repeat_proportion'] = 1.0 * df_peaks['repeat_count'] / df_peaks['length']
    df_peaks['GC_content'] = [GC(seq_record.seq) for seq_record in df_peaks['seq_records']]

    # motif features
    for fn in next(os.walk(mast_dir))[2]:
        mast_fp = os.path.join(mast_dir, fn)
        site_name = os.path.splitext(fn)[0]
        print 'adding site %s' % site_name
        df_mast_cols = _get_df_mast(mast_fp, site_name, df_peaks.drop(['sample_count_distinct', 'seq_records', 'repeat_count'], axis=1))
        df_peaks = df_peaks.join(df_mast_cols, on=('chr', 'start', 'end'))
        df_peaks.fillna(0, inplace=True)

    # background features
    k_values = (2, 3, 6)
    df_kmers = pd.DataFrame([get_kmer_dict(seq_record, k_values) for seq_record in df_peaks['seq_records']])
    df_peaks = df_peaks.join(df_kmers)
    # drop unnecessary columns
    df_peaks.drop(['sample_count_distinct', 'seq_records', 'repeat_count'], axis=1, inplace=True)

    # drop chr,start,end columns
    df_peaks.drop(['chr', 'start', 'end'], axis=1, inplace=True)
    df_peaks.to_csv(output_fp, sep='\t', index=False)


def _get_df_mast(mast_fp, site_name, df_peaks):
    """Return dataframe with aggregated P53match score and count columns.

    Join with intervals from `df_peaks`, use df_peaks interval identifier as index.

    column : aggregation function
    - P53match_score_{site_name} : max
    - P53match_count_{site_name} : len
    """
    with open(mast_fp) as f:
        comments, names = itertools.takewhile(lambda line: line.startswith('#'), f)
        names = names[1:].split()

    colname_score = 'P53match_score_' + site_name
    colname_count = 'P53match_count_' + site_name

    df_mast = pd.read_table(mast_fp, comment='#', sep='\s+', index_col=False, header=None, names=names)
    df_mast[colname_score] = -np.log10(df_mast['hit_p-value'])
    df_mast.rename(columns={'sequence_name': 'chr'}, inplace=True)

    # merge with chr, start, end from df_peaks
    df_merge_mast_all = df_peaks.merge(df_mast, on='chr')
    intersect = (df_merge_mast_all.hit_start >= df_merge_mast_all.start) & (df_merge_mast_all.hit_end <= df_merge_mast_all.end)
    df_merge_mast_intersect = df_merge_mast_all[intersect]
    group_merge_mast_intersect = df_merge_mast_intersect.groupby(list(df_peaks.columns))

    df_mast_cols = group_merge_mast_intersect.agg({'hit_p-value': len, colname_score: max}).reset_index()
    df_mast_cols[colname_count] = df_mast_cols['hit_p-value'].astype(int)
    df_mast_cols.set_index(['chr', 'start', 'end'], inplace=True)
    df_mast_cols = df_mast_cols.loc[:, (colname_score, colname_count)]
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


def get_reverse_complement(seq_str):
    """Return reverse complement of a sequence string."""
    return ''.join([reverse_complements[base] for base in seq_str[::-1]])


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged_bed', required=True,
                        help='BED file with columns chr/start/end/sample_count_distinct/max_MACS')
    parser.add_argument('--merged_fasta_softmask', required=True,
                        help='Soft-masked FASTA file for intervals in merged_bed')
    parser.add_argument('--mast_dir', required=True,
                        help='Directory containing MAST files for motifs')
    parser.add_argument('-o', required=True,
                        help='Output filepath')
    args = parser.parse_args()
    output_features(args.merged_bed, args.merged_fasta_softmask, args.mast_dir, args.o)
