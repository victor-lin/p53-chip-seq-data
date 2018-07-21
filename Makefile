bed_dir := data/SampleBEDs
fe_dir := data/FEfiles
mast_dir := data/MAST

dm6_dir := /ufrc/zhou/share/genomes/dm6
dm6_genome_fasta := $(dm6_dir)/Sequence/WholeGenomeFasta/genome.fa
dm6_genome_phastcons := $(dm6_dir)/phastCons27way/dm6.27way.phastCons.bed

MINSAMPLES := 2
REP_THRESHOLD_TYPE := max
REP_CUTOFF := 0.1
N_NONBINDING_INTERVALS := 5

TEST_SIZE := 0.5
datafile_appendix := minsamples_$(MINSAMPLES)__rep_$(REP_THRESHOLD_TYPE)$(REP_CUTOFF)__nonbinding_$(N_NONBINDING_INTERVALS)
train_out := etc/peaks_merged_features_train__$(datafile_appendix).txt
test_out := etc/peaks_merged_features_test__$(datafile_appendix).txt

etc/peaks_binding_all_samples.txt: scripts/concat_sample_beds.R $(bed_dir)/*.bed $(fe_dir)/*.xls
	Rscript $< --bed_directory $(bed_dir) --fe_directory $(fe_dir) --ignore_chr chrM -o $@

etc/peaks_binding_merged_maxMACS.bed: etc/peaks_binding_all_samples.txt
	tail -n+2 $< | sort -k1,1 -k2,2n | bedtools merge -i - -c 5,6 -o count_distinct,max > $@

etc/peaks_binding_merged_subset.txt: scripts/subset_binding_peaks.py etc/peaks_binding_merged_maxMACS.bed
	python $< --merged_bed etc/peaks_binding_merged_maxMACS.bed \
			  --minsamples $(MINSAMPLES) \
			  --rep_threshold_type $(REP_THRESHOLD_TYPE) \
			  --rep_cutoff $(REP_CUTOFF) \
			  --genome_fasta $(dm6_genome_fasta) > $@

etc/peaks_nonbinding.txt: scripts/generate_nonbinding_peaks.py etc/peaks_binding_merged_subset.txt
	python $< --binding_peaks etc/peaks_binding_merged_subset.txt \
			  --genome_fasta $(dm6_genome_fasta) \
			  --rep_threshold_type $(REP_THRESHOLD_TYPE) \
			  --rep_cutoff $(REP_CUTOFF) \
			  --num_intervals_per_side $(N_NONBINDING_INTERVALS) \
			  -o $@

etc/peaks_all.txt: etc/peaks_binding_merged_subset.txt etc/peaks_nonbinding.txt
	(head -n1 etc/peaks_binding_merged_subset.txt; \
	 ((tail -n+2 etc/peaks_binding_merged_subset.txt; \
	   tail -n+2 etc/peaks_nonbinding.txt) | sort -k1,1 -k2,2n)) > $@

etc/peaks_all.bed: etc/peaks_all.txt
	cut -f1,2,3 etc/peaks_all.txt | tail -n+2 > $@

etc/peaks_all_phastcons.bed: etc/peaks_all.bed $(dm6_genome_phastcons)
	bedmap --echo --delim "\t" --mean etc/peaks_all.bed $(dm6_genome_phastcons) > $@

etc/peaks_merged_features.txt: scripts/peak_features.py etc/peaks_all.txt etc/peaks_all_phastcons.bed $(mast_dir) $(dm6_genome_fasta)
	python $< --peaks_all etc/peaks_all.txt \
		      --peaks_all_phastcons etc/peaks_all_phastcons.bed \
			  --mast_dir $(mast_dir) \
			  --genome_fasta $(dm6_genome_fasta) \
			  -o $@

train_test_split: scripts/peak_features_preprocess.py etc/peaks_merged_features.txt
	python $< --data_file etc/peaks_merged_features.txt \
			  --test_size $(TEST_SIZE) \
			  --train_out $(train_out) \
			  --test_out $(test_out)
