import argparse
import os
import itertools
import pandas as pd
import numpy as np
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqUtils import GC


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

    # background features
    df_kmers = pd.DataFrame([get_kmer_dict(seq_record) for seq_record in df_peaks['seq_records']])
    df_peaks = df_peaks.join(df_kmers)
    # drop unnecessary columns
    df_peaks.drop(['sample_count_distinct', 'seq_records', 'repeat_count'], axis=1, inplace=True)

    # motif features
    for fn in next(os.walk(mast_dir))[2]:
        mast_fp = os.path.join(mast_dir, fn)
        site_name = os.path.splitext(fn)[0]
        print 'adding site %s' % site_name
        df_mast_cols = _get_df_mast(mast_fp, site_name, df_peaks)
        df_peaks = df_peaks.join(df_mast_cols, on=('chr', 'start', 'end'))
        df_peaks.fillna(0, inplace=True)

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


def get_kmer_dict(seq_record):
    """"Return OrderedDict of 2-mer and 3-mer counts.

    Counts are merged forward and backwards strand for reverse complement pairs.
    Remove duplicate kmers after merging reverse complements.

    seq_record : Bio.Seq.Seq object
    """
    dup_2mers = ['GG', 'GT', 'CT', 'TG', 'TC', 'TT']
    kmers_2 = [''.join(x) for x in itertools.product('GCAT', repeat=2)]
    kmers_2 = [x for x in kmers_2 if x not in dup_2mers]
    counts2 = {kmer: seq_record.seq.upper().count_overlap(kmer) + seq_record.seq.upper().reverse_complement().count_overlap(kmer)
               for kmer in kmers_2}
    sum_counts2 = sum(counts2.values())
    if sum_counts2 == 0:
        kmer_dict = OrderedDict((k, 0) for k, v in counts2.iteritems())
    else:
        kmer_dict = OrderedDict((k, 1. * v / sum_counts2) for k, v in counts2.iteritems())

    dup_3mers = ['GGG', 'GGC', 'GGT', 'GCG', 'GCT', 'GAG', 'GAT', 'GTG', 'GTC', 'GTT', 'CGG', 'CGT', 'CCT', 'CAT', 'CTG', 'CTT', 'AGT', 'ATT', 'TGG', 'TGC', 'TGA', 'TGT', 'TCG', 'TCC', 'TCT', 'TAG', 'TAC', 'TAT', 'TTG', 'TTC', 'TTA', 'TTT']
    kmers_3 = [''.join(x) for x in itertools.product('GCAT', repeat=3)]
    kmers_3 = [x for x in kmers_3 if x not in dup_3mers]
    counts3 = {kmer: seq_record.seq.upper().count_overlap(kmer) + seq_record.seq.upper().reverse_complement().count_overlap(kmer)
               for kmer in kmers_3}
    sum_counts3 = sum(counts3.values())
    if sum_counts3 == 0:
        kmer_dict.update(OrderedDict((k, 0) for k, v in counts3.iteritems()))
    else:
        kmer_dict.update(OrderedDict((k, 1. * v / sum_counts3) for k, v in counts3.iteritems()))
    return kmer_dict


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
