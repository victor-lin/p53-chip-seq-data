import os
import pandas as pd
from optparse import OptionParser

from master_table import anno_category


parser = OptionParser()
parser.add_option('--anno_directory',
                  help='BED files directory')
parser.add_option('-o', '--out_file',
                  help='output filepath')
parser.add_option('--ignore_chr',
                  help='ignore a chromosome from analysis')
(options, args) = parser.parse_args()

fnames = os.walk(options.anno_directory).next()[2]
for fname in fnames:
    fpath = os.path.join(options.anno_directory, fname)
    if fname.find('.anno') is -1:
        raise Exception('Invalid annotation file: {}'.format(fpath))
    sample_name = fname[:fname.find('.anno')]
    sample_df = pd.read_table(fpath)
    sample_df.columns = ['PeakID' if 'PeakID' in x
                         else x for x in sample_df.columns]
    sample_df['Sample Name'] = sample_name
    if (fname == fnames[0]):
        all_df = sample_df.copy()
    else:
        all_df = all_df.append(sample_df)
if options.ignore_chr:
    all_df = all_df[all_df.Chr != options.ignore_chr]

all_df['Annotation Category'] = all_df['Detailed Annotation'].apply(anno_category)
all_df[all_df.columns].to_csv(options.out_file, sep='\t', index=False)
