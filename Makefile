bed_dir := data/SampleBEDs
fe_dir := data/FEfiles
anno_dir := data/SampleAnnos
mast_bed := data/mast_out.bed

concat_bed := etc/MACSscore_summary_valid_fe.bed
merged_bed := etc/MACSscore_summary_valid_merged.bed
merged_anno := etc/MACSscore_summary_valid_merged.anno
concat_anno := etc/all_samples.anno
fasta_file := etc/target.fa

master_table_melted := results/ChIP_master_table_samples.txt


$(concat_bed): concat_sample_beds.R $(bed_dir)/*.bed $(fe_dir)/*.xls
	Rscript $< --bed_directory $(bed_dir) --fe_directory $(fe_dir) --ignore_chr chrM -o $@

$(merged_bed): $(concat_bed)
	tail -n+2 $< | sort -k1,1 -k2,2n | bedtools merge -i - > $@

$(merged_anno): $(merged_bed)
	annotatePeaks.pl $< dm6 > $@

$(fasta_file): $(merged_bed)
	homerTools extract $< /ufrc/zhou/share/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa -mask -fa > $@

$(concat_anno): concat_sample_annos.py $(anno_dir)/*.anno
	python $< --anno_directory $(anno_dir) --ignore_chr chrM -o $@

$(master_table_melted): master_table.py $(merged_bed) $(mast_bed) $(concat_bed) $(fasta_file) $(merged_anno)
	python $< --merged_file $(merged_bed) --samples_file $(concat_bed) --fasta_file $(fasta_file) --anno_file $(merged_anno) -o $@

results/ChIP_master_table_fe.txt: pivot_master_table_source.py $(merged_bed) $(mast_bed) $(concat_bed) $(fasta_file) $(merged_anno)
	python $< --merged_file $(merged_bed) --mast_file $(mast_bed) --samples_file $(concat_bed) --fasta_file $(fasta_file) --anno_file $(merged_anno) -o $@ --sample_col fe

results/ChIP_master_table_macs.txt: pivot_master_table_source.py $(merged_bed) $(mast_bed) $(concat_bed) $(fasta_file) $(merged_anno)
	python $< --merged_file $(merged_bed) --mast_file $(mast_bed) --samples_file $(concat_bed) --fasta_file $(fasta_file) --anno_file $(merged_anno) -o $@ --sample_col macs
