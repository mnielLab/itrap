#!/bin/bash

# Helle Povlsen
# 02.06.19

PLATFORM="ILLUMINA"
EXP="exp2_MHC_"$PLATFORM # Variable

PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
IN=$PARENT_DIR/data/$EXP/processed/longranger_out
OUT=$PARENT_DIR/data/$EXP/processed/longranger_clean

mkdir -p $OUT

for DIR in `find $IN/* -maxdepth 0 -type d`; do
   
(   echo $DIR
	BARCODE_TYPE=$(basename $DIR)
    #BARCODE_TYPE=`basename $DIR | sed 's/-.*//'`
    echo $BARCODE_TYPE
	
	# output from longranger
	LONGRANGER_OUT="$DIR/outs/barcoded.fastq.gz"

	#if [[ ! -e $LONGRANGER_OUT ]]; then
	#	gunzip "$LONGRANGER_OUT.gz"
	#fi

	# remove read2 (which was just added as a dummy because its required by longranger)
	zcat $LONGRANGER_OUT | paste - - - - | awk --posix '{ if (NR % 2 != 0) { print }}' | tr "\t" "\n" > $OUT/$BARCODE_TYPE.R1.fq

	# filter for reads where a barcode was found
	cat $OUT/$BARCODE_TYPE.R1.fq | paste - - - - | grep "BX:Z:" | tr "\t" "\n" > $OUT/$BARCODE_TYPE.R1.gems.fq

	# also remove the first 3 bases base on each seq, which is the last bases of the UMI (longranger removed 7 bases)
	cat $OUT/$BARCODE_TYPE.R1.gems.fq | awk --posix '{ if (NR % 2 == 0) { sub(/^.{3}/,"") }; print }' > $OUT/$BARCODE_TYPE.R1.gems.no_umi.fq

	sed -n '1~4s/^@/>/p;2~4p' $OUT/$BARCODE_TYPE.R1.gems.no_umi.fq > $OUT/$BARCODE_TYPE.R1.gems.no_umi.fa

) &

done

wait


