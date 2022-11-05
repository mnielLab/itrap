#!/bin/bash

IN_DIR=$(echo $1 | rev | cut -d'/' -f3- | rev) #$(dirname $INP)
SB_DIR=$(find $IN_DIR/Reports/html/*  -maxdepth 0 -type d | sed 's/.*-//')
FQ_FIL=$IN_DIR/$SB_DIR/TCR/*R1_001.fastq.gz

GEM_COUNTS=$2
UMI_PER_GEM=$3

echo "Checking fastq file"
zcat $FQ_FIL | head

# Extract putative GEMs (first 16nt) from raw fasta.
# The list should not be sorted because replicates of the same GEM tells us how many reads per GEM we have.
zcat $FQ_FIL | grep -E '^[ATCG]{16}' | cut -c-16 | awk '{print $0"-1"}' > $GEM_COUNTS
# list including UMI sequence (10nt)
# First get all unique combinations of GEM+UMI, then cut out the GEMs and count how many times each GEM occurs (number unique UMIs)
zcat $FQ_FIL | grep -E '^[ATCG]{26}' | cut -c-26 | sort | uniq | cut -c-16 | awk '{print $0"-1"}' | sort | uniq -c | awk '{print $2, $1}' > $UMI_PER_GEM