import argparse
import numpy as np
import pandas as pd
from bed import BedFile, Interval
from helper_functions import get_genome_seq_records_dict, peak_file_cols

rs = 0
np.random.seed(rs)


def generate_nonbinding_intervals(bed_df, genome_seq_records, maxrep,
                                  min_dist, max_dist, step,
                                  empty_columns, output_fp):
    """Return dataframe of nonbinding intervals given conditions

    Parameters
    ----------
    bed_df : pd.DataFrame
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
    output_fp : filepath
        filepath to write output to

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
    assert list(df.columns) == peak_file_cols
    df.to_csv(output_fp, sep='\t', index=False)


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--binding_peaks', required=True,
                        help='tab-delimited file with columns chr/start/end/id/...')
    parser.add_argument('--genome_fasta', required=True,
                        help='FASTA file for genome (one sequence per chromosome)')
    parser.add_argument('--maxrep', required=False, type=float, default=1.0,
                        help='maximum proportion of repeats for binding intervals to be included')
    parser.add_argument('-o', required=True,
                        help='Output filepath')
    args = parser.parse_args()

    df_binding = pd.read_table(args.binding_peaks)

    genome_seq_records = get_genome_seq_records_dict(args.genome_fasta)
    generate_nonbinding_intervals(df_binding, genome_seq_records,
                                  args.maxrep, 1000, 10000, step=10000,
                                  empty_columns=['sample_count_distinct', 'max_MACS_score'],
                                  output_fp=args.o)
