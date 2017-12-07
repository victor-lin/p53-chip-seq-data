mast_dir=data/MAST
out_file='etc/mast_concat.bed'

if [ -f $out_file ]; then
  rm $out_file
fi

for fpath in $mast_dir/*.txt
do
  echo "Processing $fpath"
  fname=${fpath##*/}
  sample=${fname%_mast.txt}

  awk -v sample="$sample" 'BEGIN { OFS = "\t"} /^#/ {next} {v=sample}{print $1, $3, $4, $6, sample}' $fpath >> "$out_file"
done
