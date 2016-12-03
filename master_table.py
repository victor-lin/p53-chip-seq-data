import sys
from pandas import read_table
from pandasql import sqldf
from optparse import OptionParser
from numpy import log10
from Bio import SeqIO


def anno_category(detailed_anno):
    """Return annotation category from detailed annotation.

    Assumes each detailed annotation has only one category name
    unless one is followed by a parenthesis.
    If so, return the category before the parenthesis.
    """
    repeats = {'RNA', 'Simple_repeat', 'Satellite', 'SINE',
               'Low_complexity', 'LTR', 'LINE', 'DNA'}
    non_repeats = {'TSS', 'TTS', 'exon', "5' UTR", "3' UTR",
                   'CpG', 'intron', 'Intergenic', 'non-coding'}
    other = {'Other', 'Unknown'}
    left_paren_index = detailed_anno.find('(')
    if left_paren_index > 0:
        detailed_anno = detailed_anno[:left_paren_index]
    for category in (repeats | non_repeats | other):
        if category in detailed_anno:
            return category


def get_anno_table(anno_file):
    """Return DataFrame from annotation file.

    Variables
    ---------
    key : columns to sort by
    other_cols : additional columns to include from the annotation file
    """
    key = ['Chr', 'Start', 'End']
    other_cols = ['Detailed Annotation', 'Gene Name']
    anno = read_table(anno_file)
    anno = anno[key + other_cols].sort_values(by=key)
    anno.reset_index(drop=True, inplace=True)
    return anno


def generate_master_table(options):
    """Generate master table containing all columns."""

    # aggregation functions for scores
    max_all_functions = {'_max': max,
                         '_all': lambda x: ','.join(map(str, x))}

    # get peaks
    peaks = read_table(options.merged_file,
                       header=None, names=('chr', 'start', 'end'))
    samples = read_table(options.samples_file)
    samples.rename(columns={col: 'sample_' + col
                            for col in ('chr', 'start', 'end', 'length')},
                   inplace=True)

    # table DataFrame
    tbl = peaks
    merged_samples = sqldf("""SELECT *
                                FROM peaks AS p
                                JOIN samples AS s
                                     ON (s.sample_start >= p.start
                                         AND s.sample_end <= p.end)""",
                           locals())
    grouped = merged_samples.groupby(['chr', 'start', 'end'])
    tbl['sample_count'] = grouped.count()['Sample_Name'].values
    for sample_name in set(samples.Sample_Name):
        tbl[sample_name] = 0

    for table_i, row in tbl.iterrows():
        sys.stdout.write('\rUpdating sample {} columns for '
                         'row {}'.format(options.sample_col,
                                         table_i + 1))
        sys.stdout.flush()
        key = tuple(list(row.values)[:3])
        group = grouped.get_group(key)
        for i, s in enumerate(group['Sample_Name']):
            if options.sample_col == 'macs':
                tbl.loc[table_i, s] = group['MACS'].iloc[i]
            elif options.sample_col == 'fe':
                tbl.loc[table_i, s] = group['FE'].iloc[i]
    print

    anno = get_anno_table(options.anno_file)
    tbl['detailed_annotation'] = anno['Detailed Annotation']
    tbl['annotation'] = tbl['detailed_annotation'].apply(anno_category)
    tbl['gene'] = anno['Gene Name']

    # aggregate samples by regions
    grouped_samples = grouped.agg({'MACS':
                                   max_all_functions,
                                   'FE':
                                   max_all_functions}).reset_index()
    grouped_samples.columns = map(''.join, grouped_samples.columns.values)

    # aggregate matrices by regions
    mast_out = read_table(options.mast_file)
    mast_out['P53match_score'] = -log10(mast_out['hit_p.value'])
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
    tbl['repeat_count'] = [seq_record.seq.count('N') for seq_record
                           in SeqIO.parse(options.fasta_file, 'fasta')]
    tbl['peak_length'] = [len(seq_record) for seq_record
                          in SeqIO.parse(options.fasta_file, 'fasta')]
    tbl['repeat_proportion'] = 1.0 * tbl['repeat_count'] / tbl['peak_length']
    return tbl


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
    master_table = generate_master_table(options)
    master_table.to_csv(options.out_file, sep='\t', index=False)
