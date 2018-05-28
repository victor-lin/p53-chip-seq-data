bed_dir := data/SampleBEDs
fe_dir := data/FEfiles
mast_dir := data/MAST

dm6_dir := /ufrc/zhou/share/genomes/dm6
dm6_chrom_sizes := $(dm6_dir)/dm6.chrom.sizes
dm6_genome_fasta := $(dm6_dir)/Sequence/WholeGenomeFasta/genome.fa
dm6_phastcons := $(dm6_dir)/phastCons27way/dm6.27way.phastCons.bed

concat_bed_macs := etc/peaks_all_samples.bed  # chr/start/end/length/sample_name/MACS_score/FE_score
merged_bed_3col := etc/peaks_merged_3col.bed  # chr/start/end
merged_bed_3col_with_negative := etc/peaks_merged_3col_with_negative.bed  # chr/start/end
merged_bed := etc/peaks_merged_maxMACS.bed  # chr/start/end/sample_count_distinct/max_MACS_score
merged_bed_with_negative := etc/peaks_merged_maxMACS_with_negative.bed  # chr/start/end/sample_count_distinct/max_MACS_score
merged_anno := etc/peaks_merged.anno
merged_fasta_mask := etc/peaks_merged_mask.fa
merged_fasta_softmask_negative := etc/peaks_merged_softmask_negative.fa
merged_phastcons_with_negative := etc/peaks_phastcons_with_negative.bed

master_table_melted := results/ChIP_master_table_samples.txt
macs_data_file := results/ChIP_MACS_merged_features.txt

$(concat_bed_macs): concat_sample_beds.R $(bed_dir)/*.bed $(fe_dir)/*.xls
	Rscript $< --bed_directory $(bed_dir) --fe_directory $(fe_dir) --ignore_chr chrM -o $@

$(merged_bed_3col): $(concat_bed_macs)
	tail -n+2 $< | sort -k1,1 -k2,2n | bedtools merge -i - > $@

$(merged_bed): $(concat_bed_macs)
	tail -n+2 $< | sort -k1,1 -k2,2n | bedtools merge -i - -c 5,6 -o count_distinct,max > $@

$(merged_bed_with_negative): generate_negative_intervals.py $(merged_bed)
	python $< --macs_bed $(merged_bed) --chrom_sizes $(dm6_chrom_sizes) | sort -k1,1 -k2,2n > $@

$(merged_bed_3col_with_negative): $(merged_bed_with_negative)
	cut -f 1,2,3 $< > $@

$(merged_anno): $(merged_bed)
	annotatePeaks.pl $< dm6 > $@

$(merged_fasta_mask): $(merged_bed_3col)
	homerTools extract $< $(dm6_genome_fasta) -mask -fa > $@

$(merged_fasta_softmask_negative): $(merged_bed_3col_with_negative)
	bedtools getfasta -fi $(dm6_genome_fasta) -bed $< -fo $@

$(merged_phastcons_with_negative): $(merged_bed_3col_with_negative) $(dm6_phastcons)
	bedmap --echo --delim "\t" --mean $< $(dm6_phastcons) > $@

$(master_table_melted): master_table.py $(merged_bed_3col) $(concat_bed_macs) $(merged_fasta_mask) $(merged_anno)
	python $< --merged_file $(merged_bed_3col) --samples_file $(concat_bed_macs) --fasta_file $(merged_fasta_mask) --anno_file $(merged_anno) -o $@

$(macs_data_file): macs_features.py $(merged_bed_with_negative) $(merged_fasta_softmask_negative) $(merged_phastcons_with_negative) $(mast_dir)
	python $< --merged_bed $(merged_bed_with_negative) --merged_fasta_softmask $(merged_fasta_softmask_negative) --merged_phastcons_bed $(merged_phastcons_with_negative) --mast_dir $(mast_dir) -o $@

results/ChIP_master_table_fe.txt: pivot_master_table.py $(master_table_melted)
	python $< --master_table $(master_table_melted) -o $@ --score fe

results/ChIP_master_table_macs.txt: pivot_master_table.py $(master_table_melted)
	python $< --master_table $(master_table_melted) -o $@ --score macs
