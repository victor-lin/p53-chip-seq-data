include make_variables

etc/peaks_binding_all_samples.txt: scripts/concat_sample_beds.py $(bed_dir)/*.bed $(fe_dir)/*.xls
	# un-merged ChIP-seq intervals from all samples
	# cols: chr/start/end/length/sample_name/MACS_score/FE_score
	python $< --bed_directory $(bed_dir) --fe_directory $(fe_dir) --ignore_chr chrM -o $@

include makefile_binding_dataset
include makefile_master_table

.PHONY: clean

clean:
	rm etc/*
