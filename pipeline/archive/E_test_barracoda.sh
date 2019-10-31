#!/usr/bin/bash


# Andrea Marquard
# 15.05.2019


# Get this script's parent dir, and work from there:
PARENT_DIR=`dirname $0`/../..

IN=$PARENT_DIR/tmp/longranger_in
OUT=$PARENT_DIR/tmp/test_barracoda

mkdir -p $OUT

# sequences:
TEMPLATESWITCH=TTTCTTATATGGG

PRIMER_R4=CGAGTACCATGGGCGTAAC

ANNEAL=GTGTGACCTTCCCCAAAAGGCGTAG



# for SAMPLE in CD8 MHC; do
for SAMPLE in MHC; do
# for SAMPLE in CD8; do

  f2="${PARENT_DIR}/tmp/barracoda_in/MHC.R1.gems.no_umi.no_adapters.trimmed.revcomp.with_sample_id.fq"
  
#   f="$IN/${SAMPLE}-2_S2_L001_R1_001.fastq"
#   f2="$OUT/$SAMPLE.raw.revcomp.fq"
#   
# 	$PARENT_DIR/scripts/utils/revcomp-fastq.pl $f > $f2


  # get fwd primer sequence
  PRIMER_FWD=$(grep -v ">" $PARENT_DIR/data/barcode_library/$SAMPLE.fwd-primer.fa)
  
  # output dir
  mkdir -p $OUT/$SAMPLE

  # # path to barracoda script
  # BARRACODA=/home/tuba/marquard/git/Barracoda/webservice
  BARRACODA=${PARENT_DIR}/scripts/barracoda

  # run it
  $BARRACODA/src/barracoda-1.0.pl  -v \
    -f $f2 \
    -m $PARENT_DIR/results/$SAMPLE.sample-key-map.tsv \
    -A $PARENT_DIR/results/$SAMPLE.sample_id.fasta \
    -B $PRIMER_FWD \
    -C 6 \
    -D $PARENT_DIR/results/$SAMPLE.a.fa \
    -E $ANNEAL \
    -F $PARENT_DIR/results/b.fa \
    -G 6 \
    -H "$TEMPLATESWITCH$PRIMER_R4" \
    -a $PARENT_DIR/results/$SAMPLE.annotations.xlsx \
    -o $OUT/$SAMPLE/out \
    -s $OUT/$SAMPLE/storage \
    -l $OUT/$SAMPLE/log \
    -k \
    $@

done






