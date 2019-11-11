#!/bin/bash

# Helle Povlsen
# 03.06.2019

# Run fastqc to see how much adapter contamination there is 
# (I've seen some small RNA RPI1 adapter)
#PLATFORM="ILLUMINA"
#EXP="exp2_MHC_"$PLATFORM # Variable

#PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
IN=$1 #$PARENT_DIR/data/$EXP/
OUT=$2 #$PARENT_DIR/data/$EXP/fastqc
IN_FIL=$3
OUT_FIL=$4

if [ -f "$IN_FIL" ]; then
    echo "$IN_FIL exist"
else 
    echo "$IN_FIL does not exist"
fi

mkdir -p $OUT

for file in `find $IN -type f \( -iname "*.fastq" -o -iname "*.fq" \)`; do
	
	FILENAME=$(basename ${file%.*})
	echo $FILENAME
	fastqc --extract -o $OUT $file
	#csplit --prefix=$OUT/$FILENAME/fastqc_data $OUT/$FILENAME/fastqc_data.txt "/>>END_MODULE/+1" "{*}"


	# Overrepresented sequences can be found in $OUT/fastqc_data08

done

wait

if [ -f "$OUT_FIL" ]; then
    echo "$OUT_FIL exist"
else 
    echo "$OUT_FIL does not exist"
fi
