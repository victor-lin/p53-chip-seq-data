import sys
from pandasql import sqldf
from optparse import OptionParser

from master_table import get_merged_peak_df, get_sample_peak_df, get_mast_data


def generate_master_table_pivoted(options):
    """Generate pivoted master table."""

    # aggregation functions for scores
    max_all_functions = {'_max': max,
                         '_all': lambda x: ','.join(map(str, x))}

    # get peaks
    peaks = get_merged_peak_df(options)
    samples = get_sample_peak_df(options.samples_file)

    # table DataFrame
    tbl = peaks
    merged_samples = sqldf("""SELECT *
                                FROM peaks AS p
                                JOIN samples AS s
                                     ON (s.sample_start >= p.start
                                         AND s.sample_end <= p.end)""",
                           locals())
    grouped = merged_samples.groupby(['chr', 'start', 'end'])
    tbl['sample_count'] = grouped.count()['sample_name'].values
    for sample_name in set(samples['sample_name']):
        tbl[sample_name] = 0

    for table_i, row in tbl.iterrows():
        sys.stdout.write('\rUpdating sample {} columns for '
                         'row {}'.format(options.sample_col,
                                         table_i + 1))
        sys.stdout.flush()
        key = tuple(list(row.values)[:3])
        group = grouped.get_group(key)
        for i, s in enumerate(group['sample_name']):
            if options.sample_col == 'macs':
                tbl.loc[table_i, s] = group['MACS'].iloc[i]
            elif options.sample_col == 'fe':
                tbl.loc[table_i, s] = group['FE'].iloc[i]
    print

    # aggregate samples by regions
    grouped_samples = grouped.agg({'MACS':
                                   max_all_functions,
                                   'FE':
                                   max_all_functions}).reset_index()
    grouped_samples.columns = map(''.join, grouped_samples.columns.values)

    # aggregate matrices by regions
    mast_out = get_mast_data(options.mast_file)
    merged_mast = sqldf("""SELECT *, mo.ROWID AS matrix_id
                             FROM peaks AS p
                             JOIN mast_out AS mo
                                  ON (mo.hit_start >= p.start
                                      AND mo.hit_end <= p.end)""",
                        locals())
    grouped2 = merged_mast.groupby(['chr', 'start', 'end'])
    # TODO: aggregate matrix_id
    grouped_mast = grouped2.agg({'P53match_score':
                                 max_all_functions}).reset_index()
    grouped_mast.columns = map(''.join, grouped_mast.columns.values)

    tbl = sqldf("""   SELECT t.*,
                             gm.P53match_score_max,
                             gs.MACS_max,
                             gs.FE_max,
                             gm.P53match_score_all,
                             gs.MACS_all,
                             gs.FE_all
                        FROM tbl AS t
                   LEFT JOIN grouped_mast AS gm
                             ON (gm.chr = t.chr
                                 AND gm.start = t.start
                                 AND gm.end = t.end)
                   LEFT JOIN grouped_samples AS gs
                             ON (gs.chr = t.chr
                                 AND gs.start = t.start
                                 AND gs.end = t.end)""",
                locals())
    tbl.to_csv(options.out_file, sep='\t', index=False)


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('--merged_file',
                      help='BED file defining merged regions')
    parser.add_option('--samples_file',
                      help='concatenated samples BED file')
    parser.add_option('--mast_file',
                      help='matrix results file')
    parser.add_option('--anno_file',
                      help='annotation file for merged regions')
    parser.add_option('--fasta_file',
                      help='FASTA file for merged regions')
    parser.add_option('-o', '--out_file',
                      help='output filepath')
    parser.add_option('--sample_col', choices=['macs', 'fe'],
                      help='value for sample column')
    (options, args) = parser.parse_args()
    generate_master_table_pivoted(options)
