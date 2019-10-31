#!/bin/bash

# Andrea Marquard
# 15.05.2019

# Get this script's parent dir, and work from there:
PARENT_DIR=`dirname $0`/../..

IN=$PARENT_DIR/tmp/longranger_out

# where to put output from this script
OUT="$PARENT_DIR/tmp/longranger_clean"
mkdir -p $OUT


for DIR in `find $IN/* -maxdepth 0 -type d`; do
   
(   echo $DIR
    SAMPLE=`basename $DIR | sed 's/-.*//'`
    echo $SAMPLE
	
	# output from longranger
	LONGRANGER_OUT="$DIR/outs/barcoded.fastq"

	if [[ ! -e $LONGRANGER_OUT ]]; then
		gunzip "$LONGRANGER_OUT.gz"
	fi

	# remove read2 (which was just added as a dummy because its required by longranger)
	cat $LONGRANGER_OUT | paste - - - - | awk --posix '{ if (NR % 2 != 0) { print }}' | tr "\t" "\n" > $OUT/$SAMPLE.R1.fq

	# filter for reads where a barcode was found
	cat $OUT/$SAMPLE.R1.fq | paste - - - - | grep "BX:Z:" | tr "\t" "\n" > $OUT/$SAMPLE.R1.gems.fq


	# also remove the first 3 bases base on each seq, which is the last bases of the UMI (longranger removed 7 bases)
	cat $OUT/$SAMPLE.R1.gems.fq | awk --posix '{ if (NR % 2 == 0) { sub(/^.{3}/,"") }; print }' > $OUT/$SAMPLE.R1.gems.no_umi.fq

) &

done

wait

# this gave 4.8 M reads in IONTORRENT.R1.gems.no_umi.fq
# of which 4.3 M contain the perfect TSO

