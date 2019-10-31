#!/bin/bash

# Get this script's parent dir, and work from there:
PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
EXP="exp3_TCR" # Variable
IN=$PARENT_DIR/data/$EXP/processed/mkfastq_in
OUT=$PARENT_DIR/data/$EXP/processed/mkfastq_out
DATASET="TCR"
INDEX="SI-GA-A2" # Variable

mkdir -p $IN
if [ ! -d "$OUT" ]; then
	mkdir -p $OUT
else
	rm -rf $OUT
	mkdir -p $OUT
fi

if [ -d "$PARENT_DIR/data/$EXP/processed/$DATASET" ]; then
	rm -rf $PARENT_DIR/data/$EXP/processed/$DATASET
fi

echo Lane,Sample,Index > $IN/index.csv
echo 1,TCR,$INDEX >> $IN/index.csv


#CELLRANGER=/home/tuba/marquard/tcr_seq/10x/apps/cellranger-3.0.0/cellranger


cellranger mkfastq --run=$PARENT_DIR/data/$EXP/raw \
                    --id=$DATASET \
                    --csv=$IN/index.csv \
                    --qc \
                    --output-dir=$OUT

mv $DATASET/ $PARENT_DIR/data/$EXP/processed/.
mv __$DATASET.* $PARENT_DIR/data/$EXP/processed/.
