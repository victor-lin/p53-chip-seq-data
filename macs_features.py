import argparse
import pandas as pd
import itertools
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

    df_peaks.drop(['chr', 'start', 'end', 'sample_count_distinct', 'seq_records', 'repeat_count'], axis=1, inplace=True)
    df_peaks.to_csv(output_fp, sep='\t', index=False)


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
    kmer_dict = OrderedDict((k, 1. * v / sum_counts2) for k, v in counts2.iteritems())

    dup_3mers = ['GGG', 'GGC', 'GGT', 'GCG', 'GCT', 'GAG', 'GAT', 'GTG', 'GTC', 'GTT', 'CGG', 'CGT', 'CCT', 'CAT', 'CTG', 'CTT', 'AGT', 'ATT', 'TGG', 'TGC', 'TGA', 'TGT', 'TCG', 'TCC', 'TCT', 'TAG', 'TAC', 'TAT', 'TTG', 'TTC', 'TTA', 'TTT']
    kmers_3 = [''.join(x) for x in itertools.product('GCAT', repeat=3)]
    kmers_3 = [x for x in kmers_3 if x not in dup_3mers]
    counts3 = {kmer: seq_record.seq.upper().count_overlap(kmer) + seq_record.seq.upper().reverse_complement().count_overlap(kmer)
               for kmer in kmers_3}
    sum_counts3 = sum(counts3.values())
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
