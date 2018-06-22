bed_dir := data/SampleBEDs
fe_dir := data/FEfiles
mast_dir := data/MAST

dm6_dir := /ufrc/zhou/share/genomes/dm6
dm6_genome_fasta := $(dm6_dir)/Sequence/WholeGenomeFasta/genome.fa
dm6_genome_phastcons := $(dm6_dir)/phastCons27way/dm6.27way.phastCons.bed

concat_bed_macs := etc/peaks_all_samples.bed  # chr/start/end/length/sample_name/MACS_score/FE_score
merged_bed := etc/peaks_merged_maxMACS.bed  # chr/start/end/sample_count_distinct/max_MACS_score
data_file_minsamples2_maxrep10 := etc/peaks_merged_features_minsamples2_maxrep10_unprocessed.txt

$(concat_bed_macs): concat_sample_beds.R $(bed_dir)/*.bed $(fe_dir)/*.xls
	Rscript $< --bed_directory $(bed_dir) --fe_directory $(fe_dir) --ignore_chr chrM -o $@

$(merged_bed): $(concat_bed_macs)
	tail -n+2 $< | sort -k1,1 -k2,2n | bedtools merge -i - -c 5,6 -o count_distinct,max > $@

$(data_file_minsamples2_maxrep10): macs_features.py $(merged_bed) $(mast_dir) $(dm6_genome_fasta) $(dm6_genome_phastcons)
	python $< --merged_bed $(merged_bed) \
			  --mast_dir $(mast_dir) \
			  --genome_fasta $(dm6_genome_fasta) \
			  --genome_phastcons $(dm6_genome_phastcons) \
			  --minsamples 2 --maxrep 0.1 \
			  -o $@
