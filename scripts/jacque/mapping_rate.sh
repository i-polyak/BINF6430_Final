#!/bin/bash

cd /scratch/caldwell.ja/final_project/kallisto

#   for loop below-- uses the files in quant/*/run_info.json 
#    d= sample number calculated above.  s = basename dir $d
#
#   jp program that sorts json files for you 
#	-- takes in directory and file names and sorts by column. 
#   after loop sort orders the rows by col4 (mapping rate)
#
# sed then creates the headder to the new file after sorting
# col just makes the whitespace separated fields into columns (padding as needed
)l


for d in quant/*/; do
  s=$(basename $d)
  jq -r --arg s "$s" '[$s, .n_processed, .n_pseudoaligned, .p_pseudoaligned] | @
tsv' ${d}run_info.json
done | sort -k4,4n > qc_mapping.tsv

sed -i '1i sample\tn_processed\tn_pseudoaligned\tp_pseudoaligned' qc_mapping.tsv
column -t qc_mapping.tsv | head -5      # lowest mapping rates
