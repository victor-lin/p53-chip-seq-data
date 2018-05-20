bed_dir := data/SampleBEDs
fe_dir := data/FEfiles
mast_dir := data/MAST

concat_bed := etc/peaks_all_samples.bed  # chr/start/end/length/sample_name/MACS_score/FE_score
merged_bed_3col := etc/peaks_merged_3col.bed  # chr/start/end
merged_bed := etc/peaks_merged_maxMACS.bed  # chr/start/end/sample_count_distinct/max_MACS_score
merged_anno := etc/peaks_merged.anno
merged_fasta_mask := etc/peaks_merged_mask.fa
merged_fasta_softmask := etc/peaks_merged_softmask.fa

master_table_melted := results/ChIP_master_table_samples.txt
master_table_new := results/ChIP_MACS_merged_features.txt

$(concat_bed): concat_sample_beds.R $(bed_dir)/*.bed $(fe_dir)/*.xls
	Rscript $< --bed_directory $(bed_dir) --fe_directory $(fe_dir) --ignore_chr chrM -o $@

$(merged_bed_3col): $(concat_bed)
	tail -n+2 $< | sort -k1,1 -k2,2n | bedtools merge -i - > $@

$(merged_bed): $(concat_bed)
	tail -n+2 $< | sort -k1,1 -k2,2n | bedtools merge -i - -c 5,6 -o count_distinct,max > $@

$(merged_anno): $(merged_bed)
	annotatePeaks.pl $< dm6 > $@

$(merged_fasta_mask): $(merged_bed_3col)
	homerTools extract $< /ufrc/zhou/share/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa -mask -fa > $@

$(merged_fasta_softmask): $(merged_bed)
	bedtools getfasta -fi /ufrc/zhou/share/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa -bed $(merged_bed) -fo $@

$(master_table_melted): master_table.py $(merged_bed_3col) $(concat_bed) $(merged_fasta_mask) $(merged_anno)
	python $< --merged_file $(merged_bed_3col) --samples_file $(concat_bed) --fasta_file $(merged_fasta_mask) --anno_file $(merged_anno) -o $@

$(master_table_new): macs_features.py $(merged_bed) $(merged_fasta_softmask) $(mast_dir)
	python $< --merged_bed $(merged_bed) --merged_fasta_softmask $(merged_fasta_softmask) --mast_dir $(mast_dir) -o $@

results/ChIP_master_table_fe.txt: pivot_master_table.py $(master_table_melted)
	python $< --master_table $(master_table_melted) -o $@ --score fe

results/ChIP_master_table_macs.txt: pivot_master_table.py $(master_table_melted)
	python $< --master_table $(master_table_melted) -o $@ --score macs