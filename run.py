#!/usr/bin/env python

import argparse
import subprocess
import os
import pandas as pd


parser = argparse.ArgumentParser(description='Create P53 figures')
parser.add_argument('-v', '--valid',
                    help='filter for valid samples')

# base dir
base_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(base_dir, 'data')
etc_dir = os.path.join(base_dir, 'etc')
result_dir = os.path.join(base_dir, 'results')

# data/
sample_beds_dir = os.path.join(data_dir, 'ValidSamples')
mast_out = os.path.join(data_dir, 'mast_out.bed')
target_fasta = os.path.join(data_dir, 'target.fa')

# etc/
samples_file = os.path.join(etc_dir, 'MACSscore_summary_valid_fe.bed')
merged_file = os.path.join(etc_dir, 'MACSscore_summary_valid_merged.bed')
anno_file = os.path.join(etc_dir, 'MACSscore_summary_valid_merged.anno')

# results/
master_table_macs_maxmat = os.path.join(result_dir, 'ChIP_master_table_macs_maxmat.txt')
master_table_fe_maxmat = os.path.join(result_dir, 'ChIP_master_table_fe_maxmat.txt')

r_cmd = 'Rscript'
r_concat_sample_beds = 'concat_sample_beds.R'
r_data_plots = 'data_plots.R'
py_master_table = 'master_table.py'

subprocess.call([r_cmd, r_concat_sample_beds,
                 '--bed_directory', sample_beds_dir,
                 '-o', samples_file,
                 '--fe_directory', os.path.join(data_dir, 'FEfiles'),
                 '--ignore', 'chrM'])
subprocess.check_call('bedtools merge -i <(tail -n+2 {} | '
                      'sort -k1,1 -k2,2n) > {}'.format(samples_file,
                                                       merged_file),
                      shell=True, executable='/bin/bash')

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
    generate_master_table(master_table_macs_maxmat, 'macs')
    generate_master_table(master_table_fe_maxmat, 'fe')

    mt_macs = pd.read_table(master_table_macs_maxmat)
    mt_fe = pd.read_table(master_table_fe_maxmat)

    mt_macs[mt_macs["Sample Count"] > 1].to_csv('ChIP_master_table_intersect_macs.txt',
                                                sep='\t', index=False)
    mt_fe[mt_fe["Sample Count"] > 1].to_csv('ChIP_master_table_intersect_fe.txt',
                                            sep='\t', index=False)
    mt_macs[mt_macs["Max Matrix Score"] > 6].to_csv('ChIP_master_table_macs_6mat.txt',
                                                    sep='\t', index=False)
    mt_fe[mt_fe["Max Matrix Score"] > 6].to_csv('ChIP_master_table_fe_6mat.txt',
                                                sep='\t', index=False)
