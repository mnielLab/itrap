#!/bin/bash

# Get this script's parent dir, and work from there:
DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
EXP="exp3_TCR" # Variable
INP=$DIR/data/$EXP/processed/mkfastq_out
OUT=$DIR/data/$EXP/processed/cellranger_out

DATASET="TCR" # Variable
PRJ=$(find $INP/Reports/html/*  -maxdepth 0 -type d | sed 's/.*-//')
REF=/home/tuba/marquard/tcr_seq/10x/apps/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0

if [ ! -d "$OUT" ]; then
        mkdir -p $OUT
else
        rm -rf $OUT
        mkdir -p $OUT
fi


cellranger vdj --id ${DATASET}_VDJ --fastqs $INP/$PRJ/$DATASET --reference $REF --sample $DATASET --chain TR

mv ${DATASET}_VDJ $OUT/

