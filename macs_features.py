import argparse


def output_features(merged_bed, merged_fasta_mask, merged_fasta_softmask, mast_dir):
    """Output data file

    Parameters
    ----------
    merged_bed : str
        filepath of BED file with columns chr/start/end/sample_count_distinct/max_MACS
    merged_fasta_mask : str
        filepath of masked FASTA file for intervals in merged_bed
    merged_fasta_softmask : str
        filepath of soft-masked FASTA file for intervals in merged_bed
    mast_dir : str
        path of directory containing MAST files for motifs
    """
    print "hello", merged_bed, merged_fasta_mask, merged_fasta_softmask, mast_dir


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged_bed', required=True,
                        help='BED file with columns chr/start/end/sample_count_distinct/max_MACS')
    parser.add_argument('--merged_fasta_mask', required=True,
                        help='Masked FASTA file for intervals in merged_bed')
    parser.add_argument('--merged_fasta_softmask', required=True,
                        help='Soft-masked FASTA file for intervals in merged_bed')
    parser.add_argument('--mast_dir', required=True,
                        help='Directory containing MAST files for motifs')
    args = parser.parse_args()
    output_features(args.merged_bed, args.merged_fasta_mask,
                    args.merged_fasta_softmask, args.mast_dir)
