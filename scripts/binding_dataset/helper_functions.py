from Bio import SeqIO

# used in subset_binding_peaks.py and  generate_nonbinding_peaks.py
peak_file_cols = ['chr', 'start', 'end', 'id', 'sample_count_distinct', 'max_MACS_score']


def get_genome_seq_records_dict(genome_fasta_file):
    return {x.name: x for x in SeqIO.parse(genome_fasta_file, 'fasta')}
