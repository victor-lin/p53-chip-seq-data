import pandas as pd
from pandasql import sqldf
from optparse import OptionParser
from numpy import log10
from Bio import SeqIO

from chip_seq import SeqSample


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
    anno = pd.read_table(anno_file)
    anno = anno[key + other_cols].sort_values(by=key)
    anno.reset_index(drop=True, inplace=True)
    return anno


def get_mast_data(mast_file):
    """Return DataFrame from MAST file."""
    mast_out = pd.read_table(mast_file)
    mast_out['P53match_score'] = -log10(mast_out['hit_p.value'])
    return mast_out


def generate_master_table_melted(options):
    """Generate master table."""
    seq_sample_attrs = ['cell_type', 'treatment_type',
                        'treatment_time', 'treatment_repeat']
    peaks = pd.read_table(options.merged_file,
                          header=None, names=('chr', 'start', 'end'))
    peaks['repeat_count'] = [seq_record.seq.count('N') for seq_record
                             in SeqIO.parse(options.fasta_file, 'fasta')]
    peaks['peak_length'] = [len(seq_record) for seq_record
                            in SeqIO.parse(options.fasta_file, 'fasta')]
    peaks['repeat_proportion'] = 1.0 * peaks['repeat_count'] / peaks['peak_length']
    samples = pd.read_table(options.samples_file)
    sample_names = set(samples['sample_name'])
    sample_obj_map = {sample_name: SeqSample(sample_name)
                      for sample_name in sample_names}
    samples.rename(columns={col: 'sample_' + col
                            for col in ('chr', 'start', 'end', 'length')},
                   inplace=True)
    tbl = sqldf("""SELECT *
                     FROM peaks AS p
                     JOIN samples AS s
                          ON (s.sample_start >= p.start
                              AND s.sample_end <= p.end)""",
                locals())
    tbl = tbl.merge(tbl['sample_name'].apply(lambda s:
                    pd.Series({attr: getattr(sample_obj_map[s], attr)
                               for attr in seq_sample_attrs})),
                    left_index=True, right_index=True)
    # TODO
    # annotation file
    # anno = get_anno_table(options.anno_file)
    # tbl['detailed_annotation'] = anno['Detailed Annotation']
    # tbl['annotation'] = tbl['detailed_annotation'].apply(anno_category)
    # tbl['gene'] = anno['Gene Name']
    # output file
    out_columns = (['chr', 'start', 'end', 'sample_name'] + seq_sample_attrs +
                   ['repeat_count', 'peak_length', 'repeat_proportion'] +
                   ['MACS_score', 'FE'])
    # + ['annotation', 'detailed_annotation', 'gene']
    tbl.to_csv(options.out_file, sep='\t', index=False,
               columns=out_columns)


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
    (options, args) = parser.parse_args()
    generate_master_table_melted(options)
