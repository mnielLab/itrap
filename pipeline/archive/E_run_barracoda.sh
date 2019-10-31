#!/bin/bash

# Andrea Marquard
# 15.05.2019

PROJECT=IONTORRENT  


# Get this script's parent dir, and work from there:
PARENT_DIR=`dirname $0`/../..

IN="$PARENT_DIR/tmp/barracoda_in"

OUT=$PARENT_DIR/results/barracoda_out



# sequences:
TEMPLATESWITCH=TTTCTTATATGGG

PRIMER_R4=CGAGTACCATGGGCGTAAC

ANNEAL=GTGTGACCTTCCCCAAAAGGCGTAG

PRIMER_FWD=GAAGTTCCAGCCAGCGTC


# output dir
mkdir -p $OUT/$SAMPLE

# # path to barracoda script
BARRACODA=/home/tuba/marquard/git/Barracoda/webservice
# BARRACODA=${PARENT_DIR}/scripts/barracoda

"$OUT/$PROJECT${IN_SUFFIX}.no_adapters.revcomp.fq"
# run it
$BARRACODA/src/barracoda-1.0.pl  -v \
  -f $IN/$PROJECT.R1.gems.no_umi.no_adapters.revcomp.fq \
  -m $PARENT_DIR/data/barcode_library/sample-key-map.txt \
  -A $PARENT_DIR/data/barcode_library/sample.fa \
  -B $PRIMER_FWD \
  -C 6 \
  -D $PARENT_DIR/data/barcode_library/a.fa \
  -E $ANNEAL \
  -F $PARENT_DIR/data/barcode_library/b.fa \
  -G 6 \
  -H "$TEMPLATESWITCH$PRIMER_R4" \
  -a $PARENT_DIR/data/barcode_library/annotations.xlsx \
  -o $OUT/$PROJECT/out \
  -s $OUT/$PROJECT/storage \
  -l $OUT/$PROJECT/log \
  -k \
  $@

    