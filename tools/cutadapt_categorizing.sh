#!/bin/bash

# Helle Povlsen
# 03.06.2019

EXP="exp3_MHC" # Variable
PROJECT="no_adapters"

PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
IN_DIR=$PARENT_DIR/data/$EXP/processed/longranger_clean
OUT_DIR=$PARENT_DIR/data/$EXP/cutadapt_naming/$PROJECT


SUFFIX="IONTORRENT.R1.gems.no_umi."$PROJECT

JELLYFISH="/home/tuba-nobackup/shared/jellyfish-2.3.0/bin/jellyfish"

# OBS!
# Run categorize_reads_by_annotations.ipynb

# categorizing fa files
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/categorized_fa
mkdir -p $OUT_DIR/jellyfish

for i in 0 1 2 3 4
do
	echo $i

	file_A=$OUT_DIR/headers/match_${i}.lst
	file_B=$OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.fa
	cat $file_B | paste - - > $OUT_DIR/tmp.tab

	# optimize for memory usage
	nlines=$(wc -l < $file_A)
	chunk=10000
	echo $nlines $chunk
	for(( l=1; l < nlines; l += chunk)) 
	do
		sed -n $l,+$((chunk - 1))p $file_A | grep -F -f - $OUT_DIR/tmp.tab | tr "\t" "\n"
	done > $OUT_DIR/categorized_fa/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.${i}.fa

	#grep -A 1 -f $OUT_DIR/headers/match_$i.lst $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.fa > $OUT_DIR/categorized_fa/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.$i.fa

	# Jellyfish
	echo "Run jellyfish"

	$JELLYFISH count -m 22 -o $OUT_DIR/jellyfish/output.${i}.jf -c 3 -s 10000000 -t 4 $OUT_DIR/categorized_fa/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.${i}.fa
	#$JELLYFISH merge -o $OUT_DIR/jellyfish/output_$i/output.jf $OUT_DIR/jellyfish/output_$i/output\_*
	$JELLYFISH dump -c $OUT_DIR/jellyfish/output.${i}.jf | sort -rn -k2 > $OUT_DIR/jellyfish/kmer_counts.${i}.tsv


done

rm $OUT_DIR/tmp.tab



