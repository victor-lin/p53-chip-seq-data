etc/MACSscore_summary_valid_fe.bed: concat_sample_beds.R data/SampleBEDs/*.bed data/FEfiles/*.xls
	Rscript concat_sample_beds.R --bed_directory data/SampleBEDs --fe_directory data/FEfiles --ignore chrM -o etc/MACSscore_summary_valid_fe.bed

etc/MACSscore_summary_valid_merged.bed: etc/MACSscore_summary_valid_fe.bed
	bedtools merge -i <(tail -n+2 etc/MACSscore_summary_valid_fe.bed | sort -k1,1 -k2,2n) > etc/MACSscore_summary_valid_merged.bed

etc/MACSscore_summary_valid_merged.anno: etc/MACSscore_summary_valid_merged.bed
	annotatePeaks.pl <(cut -f 1,2,3 etc/MACSscore_summary_valid_merged.bed) dm6 > etc/MACSscore_summary_valid_merged.anno

etc/all_samples.anno: concat_sample_annos.py data/SampleAnnos/*.anno
	python concat_sample_annos.py --anno_directory data/SampleAnnos --ignore chrM -o etc/all_samples.anno

results/ChIP_master_table.txt: master_table.py etc/MACSscore_summary_valid_merged.bed data/mast_out.bed etc/MACSscore_summary_valid_fe.bed data/target.fa etc/MACSscore_summary_valid_merged.anno
	python master_table.py --merged_file etc/MACSscore_summary_valid_merged.bed --mast_file data/mast_out.bed --samples_file etc/MACSscore_summary_valid_fe.bed --fasta_file data/target.fa --anno_file etc/MACSscore_summary_valid_merged.anno -o results/ChIP_master_table.txt
