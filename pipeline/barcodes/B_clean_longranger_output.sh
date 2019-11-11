#!/bin/bash

# Helle Povlsen
# 02.06.19

PLATFORM=$5 #"ILLUMINA"
#EXP="exp2_MHC_"$PLATFORM # Variable

#PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
IN=$2 #$PARENT_DIR/data/$EXP/processed/longranger_out
OUT_FIL=$3 #$PARENT_DIR/data/$EXP/processed/longranger_clean
OUT=$4

mkdir -p $OUT

# output from longranger
LONGRANGER_OUT=$1 #"$IN/$PLATFORM/outs/barcoded.fastq.gz"

#if [[ ! -e $LONGRANGER_OUT ]]; then
#	gunzip "$LONGRANGER_OUT.gz"
#fi

if [ "$OUT_FIL" = "$OUT/$PLATFORM.R1.gems.no_umi.fq" ]; then
	echo $OUT_FIL
else
	echo Mismatch in filenames
	echo $OUT_FIL
	echo $OUT/$PLATFORM.R1.gems.no_umi.fq
fi


# remove read2 (which was just added as a dummy because its required by longranger)
zcat $LONGRANGER_OUT | paste - - - - | awk --posix '{ if (NR % 2 != 0) { print }}' | tr "\t" "\n" > $OUT/$PLATFORM.R1.fq

# filter for reads where a barcode was found
cat $OUT/$PLATFORM.R1.fq | paste - - - - | grep "BX:Z:" | tr "\t" "\n" > $OUT/$PLATFORM.R1.gems.fq

# also remove the first 3 bases base on each seq, which is the last bases of the UMI (longranger removed 7 bases)
cat $OUT/$PLATFORM.R1.gems.fq | awk --posix '{ if (NR % 2 == 0) { sub(/^.{3}/,"") }; print }' > $OUT/$PLATFORM.R1.gems.no_umi.fq

sed -n '1~4s/^@/>/p;2~4p' $OUT/$PLATFORM.R1.gems.no_umi.fq > $OUT/$PLATFORM.R1.gems.no_umi.fa #$OUT_FIL #



wait


