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


def get_repeat_type(detailed_anno):
    """Return annotation subsection repeat type from detailed annotation.

    Return the name after the rightmost vertical bar (if present).
    If no bar is present, then return None.
    """
    last_vbar_index = detailed_anno.rfind('|')
    if last_vbar_index > 0:
        return detailed_anno[last_vbar_index + 1:]


def get_anno_df(filepath, key=None, other_cols=None):
    """Return DataFrame from HOMER annotation output file.

    Parameters
    ----------
    filepath : filepath
    key : columns to sort by
    other_cols : additional columns to include from the annotation file
    """
    if not key:
        key = ['Chr', 'Start', 'End']
    if not other_cols:
        other_cols = ['Detailed Annotation', 'Gene Name']
    anno = pd.read_table(filepath)
    anno = anno[key + other_cols].sort_values(by=key)
    anno.reset_index(drop=True, inplace=True)
    return anno


def get_mast_data(mast_file):
    """Return DataFrame from MAST file.

    Generate P53match_score from hit_p.value.
    """
    mast_out = pd.read_table(mast_file)
    mast_out['P53match_score'] = -log10(mast_out['hit_p.value'])
    return mast_out


def get_merged_peak_df(options):
    """Return pandas.DataFrame of merged peak data.

    BED columns:
    * chr
    * start
    * end

    Repeat info columns:
    * repeat_count
    * peak_length
    * repeat_proportion

    Annotation data columns:
    * detailed_annotation
    * annotation
    * gene

    """
    peaks = pd.read_table(options.merged_file,
                          header=None, names=('chr', 'start', 'end'))
    # add repeat information
    peaks['repeat_count'] = [seq_record.seq.count('N') for seq_record
                             in SeqIO.parse(options.fasta_file, 'fasta')]
    peaks['peak_length'] = [len(seq_record) for seq_record
                            in SeqIO.parse(options.fasta_file, 'fasta')]
    peaks['repeat_proportion'] = 1.0 * peaks['repeat_count'] / peaks['peak_length']
    # annotation data
    anno = get_anno_df(options.anno_file)
    peaks['detailed_annotation'] = anno['Detailed Annotation']
    peaks['annotation'] = peaks['detailed_annotation'].apply(anno_category)
    peaks['repeat_type'] = peaks['detailed_annotation'].apply(get_repeat_type)
    peaks['gene'] = anno['Gene Name']
    return peaks


def get_sample_peak_df(options, seq_sample_attrs):
    """Return pandas.DataFrame of sample peak data.

    Add additional columns for seq_sample_attrs.
    """
    samples = pd.read_table(options.samples_file)
    sample_names = set(samples['sample_name'])
    sample_obj_map = {sample_name: SeqSample(sample_name)
                      for sample_name in sample_names}
    samples.rename(columns={col: 'sample_' + col
                            for col in ('chr', 'start', 'end', 'length')},
                   inplace=True)
    samples = samples.merge(samples['sample_name'].apply(lambda s:
                            pd.Series({attr: getattr(sample_obj_map[s], attr)
                                       for attr in seq_sample_attrs})),
                            left_index=True, right_index=True)
    return samples


def generate_master_table_melted(options):
    """Generate master table."""
    seq_sample_attrs = ['cell_type', 'treatment_type',
                        'treatment_time', 'treatment_repeat']
    peaks = get_merged_peak_df(options)
    # get sample info
    samples = get_sample_peak_df(options, seq_sample_attrs)
    tbl = sqldf("""SELECT *
                     FROM peaks AS p
                     JOIN samples AS s
                          ON (s.sample_start >= p.start
                              AND s.sample_end <= p.end)""",
                locals())
    # TODO: add constraint s.sample_chr = p.chr ?
    # TODO: add MAST file data
    # output file
    out_columns = (['chr', 'start', 'end', 'sample_name'] + seq_sample_attrs +
                   ['repeat_count', 'peak_length', 'repeat_proportion'] +
                   ['MACS_score', 'FE'] +
                   ['detailed_annotation', 'annotation', 'repeat_type', 'gene'])
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
