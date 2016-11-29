#!/usr/bin/env python

import subprocess
import os
import pandas as pd


# base dir
base_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(base_dir, 'data')
etc_dir = os.path.join(base_dir, 'etc')
result_dir = os.path.join(base_dir, 'results')

# data/
sample_beds_dir = os.path.join(data_dir, 'SampleBEDs')
sample_annos_dir = os.path.join(data_dir, 'SampleAnnos')
mast_out = os.path.join(data_dir, 'mast_out.bed')
target_fasta = os.path.join(data_dir, 'target.fa')

# etc/
samples_file = os.path.join(etc_dir, 'MACSscore_summary_valid_fe.bed')
merged_file = os.path.join(etc_dir, 'MACSscore_summary_valid_merged.bed')
anno_file = os.path.join(etc_dir, 'MACSscore_summary_valid_merged.anno')
sample_annos_file = os.path.join(etc_dir, 'all_samples.anno')

# results/
mt_name = 'ChIP_master_table_{}_{}.txt'
mt_macs_maxmat = os.path.join(result_dir, mt_name.format('macs', 'maxmat'))
mt_fe_maxmat = os.path.join(result_dir, mt_name.format('fe', 'maxmat'))
mt_macs_intersect = os.path.join(result_dir, mt_name.format('macs', 'intersect'))
mt_fe_intersect = os.path.join(result_dir, mt_name.format('fe', 'intersect'))
mt_macs_6mat = os.path.join(result_dir, mt_name.format('macs', '6mat'))
mt_fe_6mat = os.path.join(result_dir, mt_name.format('fe', '6mat'))

r_cmd = 'Rscript'
r_concat_sample_beds = 'concat_sample_beds.R'
py_concat_sample_annos = 'concat_sample_annos.py'
r_data_plots = 'data_plots.R'
py_master_table = 'master_table.py'


def setup():
    subprocess.call([r_cmd, r_concat_sample_beds,
                     '--bed_directory', sample_beds_dir,
                     '-o', samples_file,
                     '--fe_directory', os.path.join(data_dir, 'FEfiles'),
                     '--ignore', 'chrM'])
    subprocess.check_call('bedtools merge -i <(tail -n+2 {} | '
                          'sort -k1,1 -k2,2n) > {}'.format(samples_file,
                                                           merged_file),
                          shell=True, executable='/bin/bash')

    subprocess.call(['python', py_concat_sample_annos,
                     '--anno_directory', sample_annos_dir,
                     '-o', sample_annos_file,
                     '--ignore', 'chrM'])

    # TODO: run annotatePeaks.pl
    # for now, assume it's in etc.
    # subprocess.call(['module', 'load', 'gcc/5.2.0', 'homer/4.8'])


def generate_master_table(out_fpath, sample_col='macs'):
    mt_args = ['--merged_file', merged_file,
               '--mast_file', mast_out,
               '--samples_file', samples_file,
               '--fasta_file', target_fasta,
               '--anno_file', anno_file,
               '--sample_col', sample_col,
               '-o', out_fpath]
    subprocess.call(['python', '-u', py_master_table] + mt_args)

if __name__ == '__main__':
    setup()
    generate_master_table(mt_macs_maxmat, 'macs')
    generate_master_table(mt_fe_maxmat, 'fe')

    mt_macs = pd.read_table(mt_macs_maxmat)
    mt_fe = pd.read_table(mt_fe_maxmat)

    mt_macs[mt_macs["sample_count"] > 1].to_csv(mt_macs_intersect,
                                                sep='\t', index=False)
    mt_fe[mt_fe["sample_count"] > 1].to_csv(mt_fe_intersect,
                                            sep='\t', index=False)
    mt_macs[mt_macs["P53match_score_max"] > 6].to_csv(mt_macs_6mat,
                                                      sep='\t', index=False)
    mt_fe[mt_fe["P53match_score_max"] > 6].to_csv(mt_fe_6mat,
                                                  sep='\t', index=False)
