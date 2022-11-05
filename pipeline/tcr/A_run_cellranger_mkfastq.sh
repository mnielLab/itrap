#!/bin/bash

# Get this script's parent dir, and work from there:
#PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
#EXP="exp6_TCR" # Variable
FQ_IN=$1 #Raw directory
#DATASET=$2 #"TCR"
#INDEX=$2 #"SI-GA-A4" # Variable
INDEX_FILE=$2 #$PARENT_DIR/data/$EXP/processed/mkfastq/index.csv

############################
# https://kb.10xgenomics.com/hc/en-us/articles/115003831603-What-to-do-if-I-forgot-which-sample-indices-I-used-during-library-construction-
# https://s3-us-west-2.amazonaws.com/10x.files/supp/cell-exp/chromium-shared-sample-indexes-plate.csv
############################

echo $FQ_IN
head $INDEX_FILE


OUT_DIR=$(echo $3 | rev | cut -d'/' -f3- | rev)
echo $OUT_DIR

SUP_DIR=$(echo $INDEX_FILE | rev | cut -d'/' -f3- | rev)
echo $SUP_DIR

cd $SUP_DIR

# If mkfastq has been run before there will exist a TCR folder that halts a rerun
# To enable rerun the TCR folder must be removed.
if [ -d "$SUP_DIR/TCR/" ]; then
        rm -rf $SUP_DIR/TCR/
fi

# Not sure what this does?
if [ -d "$OUT_DIR" ]; then
	rm -rf $OUT_DIR
	mkdir $OUT_DIR
fi

#mkdir -p $IN
#if [ ! -d "$OUT" ]; then
#	mkdir -p $OUT
#else
#	rm -rf $OUT
#	mkdir -p $OUT
#fi
#
#if [ -d "$PARENT_DIR/data/$EXP/processed/$DATASET" ]; then
#	rm -rf $PARENT_DIR/data/$EXP/processed/$DATASET
#fi

#echo Lane,Sample,Index > $INDEX_FILE #$IN/index.csv
#echo 1,TCR,$INDEX >> $INDEX_FILE #$IN/index.csv

#CELLRANGER=/home/tuba/marquard/tcr_seq/10x/apps/cellranger-3.0.0/cellranger
#CELLRANGER=/home/tuba-nobackup/shared/cellranger/cellranger-3.1.0/cellranger
CELLRANGER=/home/tuba-nobackup/shared/cellranger/cellranger-6.0.0/cellranger
$CELLRANGER mkfastq --run=$FQ_IN \
                    --id=TCR \
                    --csv=$INDEX_FILE \
                    --output-dir=$OUT_DIR


#SUP_DIR=$(echo $INDEX_FILE | rev | cut -d'/' -f3- | rev)
#

#
#mv TCR/ $SUP_DIR/.
#mv __TCR.* $SUP_DIR/.
