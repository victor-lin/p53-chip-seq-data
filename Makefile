bed_dir := data/SampleBEDs
fe_dir := data/FEfiles
anno_dir := data/SampleAnnos
mast_bed := data/mast_out.bed
fasta_file := data/target.fa

concat_bed := etc/MACSscore_summary_valid_fe.bed
merged_bed := etc/MACSscore_summary_valid_merged.bed
merged_anno := etc/MACSscore_summary_valid_merged.anno
concat_anno := etc/all_samples.anno

master_table_melted := results/ChIP_master_table.txt


$(concat_bed): concat_sample_beds.R $(bed_dir)/*.bed $(fe_dir)/*.xls
ifdef ignore
	Rscript $< --bed_directory $(bed_dir) --fe_directory $(fe_dir) --ignore_chr chrM --ignore $(ignore) -o $@
else
	Rscript $< --bed_directory $(bed_dir) --fe_directory $(fe_dir) --ignore_chr chrM -o $@
endif

$(merged_bed): $(concat_bed)
	tail -n+2 $< | sort -k1,1 -k2,2n | bedtools merge -i - > $@

$(merged_anno): $(merged_bed)
	cut -f 1,2,3 $< | annotatePeaks.pl dm6 - > $@

$(concat_anno): concat_sample_annos.py $(anno_dir)/*.anno
ifdef ignore
	python $< --anno_directory $(anno_dir) --ignore_chr chrM --ignore $(ignore) -o $@
else
	python $< --anno_directory $(anno_dir) --ignore_chr chrM -o $@
endif

$(master_table_melted): master_table.py $(merged_bed) $(mast_bed) $(concat_bed) $(fasta_file) $(merged_anno)
	python $< --merged_file $(merged_bed) --mast_file $(mast_bed) --samples_file $(concat_bed) --fasta_file $(fasta_file) --anno_file $(merged_anno) -o $@

results/ChIP_master_table_fe.txt: pivot_master_table_source.py $(merged_bed) $(mast_bed) $(concat_bed) $(fasta_file) $(merged_anno)
	python $< --merged_file $(merged_bed) --mast_file $(mast_bed) --samples_file $(concat_bed) --fasta_file $(fasta_file) --anno_file $(merged_anno) -o $@ --sample_col fe

results/ChIP_master_table_macs.txt: pivot_master_table_source.py $(merged_bed) $(mast_bed) $(concat_bed) $(fasta_file) $(merged_anno)
	python $< --merged_file $(merged_bed) --mast_file $(mast_bed) --samples_file $(concat_bed) --fasta_file $(fasta_file) --anno_file $(merged_anno) -o $@ --sample_col macs
