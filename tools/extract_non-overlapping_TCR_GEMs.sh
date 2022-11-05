#!/bin/bash

KMA_OUT=$1 #mapping.frag.gz
CEL_OUT=$2 #all_contig_annotations.csv
TCR_DIR=$3
UMI_OUT=$4
GEM_OUT=$5

IN_DIR=$(echo $TCR_DIR | rev | cut -d'/' -f3- | rev)
OUT_DIR=$(echo $UMI_OUT | rev | cut -d'/' -f2- | rev)
echo $UMI_OUT
echo $OUT_DIR

SUB_DIR=$(find $IN_DIR/Reports/html/*  -maxdepth 0 -type d | sed 's/.*-//')
TCR_FIL=$IN_DIR/$SUB_DIR/TCR/*R1_001.fastq.gz

echo $SUB_DIR
zcat $TCR_FIL | head

# Get the set of GEMs that contains BC, but were not annotated with cellranger
tail -n +2 $CEL_OUT | cut -d',' -f1 | sort | uniq > $OUT_DIR/tcr_gems.lst
zcat $KMA_OUT | cut -f7 | cut -d':' -f5 | sort | uniq > $OUT_DIR/mapped_gems.lst

echo "tcr_gems.lst"
head $OUT_DIR/tcr_gems.lst
echo "mapped_gems.lst"
head $OUT_DIR/mapped_gems.lst

grep -vf $OUT_DIR/tcr_gems.lst $OUT_DIR/mapped_gems.lst | cut -d'-' -f1 > $OUT_DIR/mapped_gems.nooverlap.seq.lst

echo "mapped_gems.nooverlap.seq.lst"
head $OUT_DIR/mapped_gems.nooverlap.seq.lst

echo "start grepping TCR reads"
zcat $TCR_FIL | grep -E '^[ATCG]{16}' > $OUT_DIR/tcr_reads.seq
echo "finalized grepping TCR reads"
head $OUT_DIR/tcr_reads.seq

echo "start grepping for BC annotated GEMs"
#grep -F -f $OUT_DIR/mapped_gems.nooverlap.seq.lst $OUT_DIR/tcr_reads.seq > $OUT_DIR/tcr_gems.nooverlap.seq
split -l 28 $OUT_DIR/mapped_gems.nooverlap.seq.lst $OUT_DIR/mapped_gems.nooverlap.seq.lst.split.
for CHUNK in $OUT_DIR/mapped_gems.nooverlap.seq.lst.split.* ; do
	idx="${CHUNK##*.}"
    grep -F -f "$CHUNK" $OUT_DIR/tcr_reads.seq > $OUT_DIR/tcr_gems.nooverlap.seq.split.$idx #& Test later if this will work!
done
# wait for all the subprocess to finish and then delete temporary files
wait
cat $OUT_DIR/tcr_gems.nooverlap.seq.split.* > $OUT_DIR/tcr_gems.nooverlap.seq
rm $OUT_DIR/mapped_gems.nooverlap.seq.lst.split.*
rm $OUT_DIR/tcr_gems.nooverlap.seq.split.*
echo "finalized"
head $OUT_DIR/tcr_gems.nooverlap.seq

# list including UMI sequence (10nt)
# First get all unique combinations of GEM+UMI, then cut out the GEMs and count how many times each GEM occurs (number unique UMIs)
cut -c-26 $OUT_DIR/tcr_gems.nooverlap.seq | sort | uniq | cut -c-16 | awk '{print $0"-1"}' | sort | uniq -c | awk '{print $2, $1}' > $UMI_OUT
# Extract putative GEMs (first 16nt) from raw fasta.
# The list should not be sorted because replicates of the same GEM tells us how many reads per GEM we have.
cut -c-16 $OUT_DIR/tcr_gems.nooverlap.seq | awk '{print $0"-1"}' | sort | uniq > $GEM_OUT
