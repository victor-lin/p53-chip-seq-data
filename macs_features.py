import argparse
import subprocess
import os
import itertools
import pandas as pd
import numpy as np
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqUtils import GC

from bed import BedFile, Interval

reverse_complements = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
genome_seq_records = None  # force this variable to be global

rs = 0  # random seed
np.random.seed(rs)


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


def generate_nonbinding_intervals(bed_df, genome_seq_records, maxrep,
                                  min_dist, max_dist, step, empty_columns):
    """Return dataframe of nonbinding intervals given conditions

    Parameters
    ----------
    bed_df : pd.Dataframe
        dataframe of binding intervals (chr/start/end)
    genome_seq_records : dict of {chr_name: Bio.Seq.Seq}
        dict of per-chromosome sequence objects
    maxrep : float
        Maximum proportion of repeats for binding intervals to be included
    min_dist : int
        minimum distance from binding interval to select from
    max_dist : int
        maximum distance from binding interval to select from
    step : int
        increment for `min_dist` and `max_dist` when no possible negative intervals are available
    empty_columns : list of str
        columns to be included with NA values

    NOTE: id is assumed to be range(0, len(bed_df))
    """
    bf = BedFile(bed_df)
    df = pd.DataFrame(columns=['chr', 'start', 'end', 'id'])

    chrom_sizes_dict = {chrom: len(genome_seq_records[chrom])
                        for chrom in genome_seq_records}

    for i, interval in enumerate(bf.intervals):
        print 'interval {} out of {}'.format(i + 1, len(bf.intervals))
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        length = interval.length

        lower = max(0, start - max_dist)
        upper = max(0, start - min_dist)
        left = None
        while not left:
            print '\tsearching in ({}, {})'.format(lower, upper)
            left = generate_random_interval(chrom, lower, upper, length, genome_seq_records,
                                            maxrep, avoid_intervals=bf.intervals)
            lower = max(0, lower - step)
            upper = max(0, upper - step)
            if lower == upper:
                print 'reached chromosome boundary. cannot search next iteration'
                break
        if left:
            df.loc[len(df)] = [left.chrom, left.start, left.end, i]

        lower = min(chrom_sizes_dict[chrom], end + min_dist)
        upper = min(chrom_sizes_dict[chrom], end + max_dist)
        right = None
        while not right:
            print '\tsearching in ({}, {})'.format(lower, upper)
            right = generate_random_interval(chrom, lower, upper, length, genome_seq_records,
                                             maxrep, avoid_intervals=bf.intervals)
            lower = max(0, lower + step)
            upper = max(0, upper + step)
            if lower == upper:
                print 'reached chromosome boundary. cannot search next iteration'
                break
        if right:
            df.loc[len(df)] = [right.chrom, right.start, right.end, i]

    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['id'] = df['id'].astype(int)
    for col_name in empty_columns:
        df[col_name] = 0
    return df


def generate_random_interval(chrom, lower, upper, length, genome_seq_records, maxrep=1.0,
                             avoid_intervals=None, min_distance=None):
    """Generate random interval given constraints.

    min_distance requires avoid_intervals
    """

    # TODO: check parameter values
    # TODO: implement min_distance

    possible_intervals = set([Interval(chrom, start, start + length)
                              for start in range(lower, upper - length)])

    possible_intervals_copy = set(possible_intervals)
    for interval in possible_intervals_copy:
        # check overlaps
        for avoid_interval in avoid_intervals:
            if interval.overlaps(avoid_interval):
                possible_intervals.remove(interval)
                break
        if interval not in possible_intervals:
            continue

        # check repeat proportion
        seq_record = genome_seq_records[chrom].seq[interval.start:interval.end]
        rep_prop = 1. * sum(map(seq_record.count, ['g', 'c', 'a', 't'])) / length
        if rep_prop > maxrep:
            possible_intervals.remove(interval)

    if not possible_intervals:
        return None
    return np.random.choice(list(possible_intervals))


def add_seq_features(df_all):
    """Requires dataframe columns chr/start/end
    Adds columns for length, repeat_proportion, GC_content
    """
    df_all['length'] = df_all['end'] - df_all['start']
    df_all['seq_records'] = [genome_seq_records[interval['chr']][interval['start']:interval['end']]
                             for i, interval in df_all.iterrows()]
    df_all['repeat_count'] = [sum(map(seq_record.seq.count, ['g', 'c', 'a', 't']))
                              for seq_record in df_all['seq_records']]
    df_all['repeat_proportion'] = 1.0 * df_all['repeat_count'] / df_all['length']
    df_all['GC_content'] = [GC(seq_record.seq) for seq_record in df_all['seq_records']]
    return df_all


def output_3col_bed(df, filepath):
    df.loc[:, ('chr', 'start', 'end')].to_csv(filepath, sep='\t', index=False, header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged_bed', required=True,
                        help='BED file with columns chr/start/end/sample_count_distinct/max_MACS')
    parser.add_argument('--mast_dir', required=True,
                        help='Directory containing MAST files for motifs')
    parser.add_argument('--genome_fasta', required=True,
                        help='FASTA file for genome (one sequence per chromosome)')
    parser.add_argument('--genome_phastcons_file', required=True,
                        help='BED file with columns chr/start/end/phastCon')
    parser.add_argument('--minsamples', required=False, type=int, default=1,
                        help='minimum number of samples for binding intervals to be included (does not apply to non-binding intervals)')
    parser.add_argument('--maxrep', required=False, type=float, default=1.0,
                        help='maximum proportion of repeats for binding intervals to be included')
    parser.add_argument('-o', required=True,
                        help='Output filepath')
    args = parser.parse_args()

    print 'reading genome FASTA...'
    genome_seq_records = {x.name: x for x in SeqIO.parse(args.genome_fasta, 'fasta')}

    print 'reading binding intervals...'
    df_binding = pd.read_table(args.merged_bed, header=None,
                               names=('chr', 'start', 'end', 'sample_count_distinct', 'max_MACS_score'))
    # add 'id' column
    df_binding.insert(3, 'id', df_binding.index)
    # add sequence features
    df_binding = add_seq_features(df_binding)
    # subset by minsamples and maxrep
    rows_binding_drop = (df_binding['sample_count_distinct'] < args.minsamples) | (df_binding['repeat_proportion'] > args.maxrep)
    df_binding = df_binding[~rows_binding_drop].copy()

    print 'generating nonbinding intervals...'
    df_nonbinding = generate_nonbinding_intervals(df_binding,
                                                  genome_seq_records,
                                                  args.maxrep, 1000, 10000, step=10000,
                                                  empty_columns=['sample_count_distinct', 'max_MACS_score'])
    # add sequence features
    df_nonbinding = add_seq_features(df_nonbinding)

    df_all = df_binding.append(df_nonbinding)

    print 'adding k-mer features'
    k_values = (2, 3)
    df_kmers = pd.DataFrame([get_kmer_dict(seq_record, k_values) for seq_record in df_all['seq_records']])
    df_all = pd.concat([df_all.reset_index(drop=True), df_kmers], axis=1)

    print 'adding phastcons...'
    bed_3col_file = 'etc/tmp_3col.bed'
    bed_3col_with_avgcons = 'etc/tmp_3col_phastcons.bed'
    output_3col_bed(df_all, bed_3col_file)
    cmd_list = ['bedmap', '--echo', '--delim "\t"', '--mean',
                bed_3col_file, args.genome_phastcons_file, '>', bed_3col_with_avgcons]
    subprocess.call(' '.join(cmd_list))

    df_phastcons = pd.read_table(bed_3col_with_avgcons, header=None,
                                 na_values='NAN', keep_default_na=False,
                                 names=('chr', 'start', 'end', 'average_phastCon'))
    assert len(df_all) == len(df_phastcons)
    df_all['average_phastCon'] = df_phastcons['average_phastCon']

    print 'adding MAST scores...'
    for fn in next(os.walk(args.mast_dir))[2]:
        mast_fp = os.path.join(args.mast_dir, fn)
        site_name = os.path.splitext(fn)[0]
        print '\tadding site %s' % site_name
        df_mast_cols = _get_df_mast(mast_fp, site_name, df_all.loc[:, ['chr', 'start', 'end']])
        df_all = df_all.join(df_mast_cols, on=('chr', 'start', 'end'))

    # drop unnecessary columns, fill NA with 0
    df_all.drop(['chr', 'start', 'end', 'seq_records', 'repeat_count'],
                axis=1, inplace=True)
    df_all.fillna(0, inplace=True)
    df_all.to_csv(args.o, sep='\t', index=False)
