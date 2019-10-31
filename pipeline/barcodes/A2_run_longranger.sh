#!/bin/bash

# Helle Povlsen
# 02.06.19

# iontorrent reads start with part of the 10X-adapter Akey, which is:
#                    CCATCTCATCCCTGCGTGTCTCCGACTCAGCACGACGCTCTTCCGATC
# The reads begin like so:                         CACGACGCTCTTCCGATCT
# Then comes 26 bases (GEM barcode, and 16bp UMI)                     NNNNNNNNNNNNNNNNNNNNNNNNNN                       
# Then the TSO (# 4.5 M reads contain the 10x TSO)                                              TTTCTTATATGGG

PLATFORM="ILLUMINA"
EXP="exp2_MHC_"$PLATFORM # Variable

PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
IN=$PARENT_DIR/data/$EXP/processed/cut_adapter
OUT=$PARENT_DIR/data/$EXP/processed/longranger_out

mkdir -p $OUT

for file in $IN/${PLATFORM}_S1_L001_R1_001.fastq ; do
	FQ1=$file
	echo $FQ1
	if [ -f "$FQ1" ]; then
		FQ2=${FQ1/_R1_/_R2_}
		echo $FQ2
	else
		echo "$FQ1 does not exist"
	fi

	# make a dummy read2 file
	cat $FQ1 | awk --posix '{ if (NR % 4 == 0) { sub(/.*/,"CCCCCCCC")} else if (NR % 2 == 0) { sub(/.*/,"NNNNNNNN")}; print}' > $FQ2

	# IMPORTANT, because if dir already exists, longranger wont run again, and we might not notice
	rm -rf $PLATFORM #$BARCODE_TYPE       # where longranger will put it
	rm -rf $OUT/$PLATFORM #$BARCODE_TYPE  # where we will put it in the end

	# Run longranger to demultplex single cells/gems
	echo "Running Longranger on $PLATFORM..." #BARCODE_TYPE
	longranger basic --id=$PLATFORM --fastqs=$IN --sample=$PLATFORM > $OUT/$PLATFORM.log 2>&1

	# Longranger places the files locally
	[[ -d $PLATFORM ]] && mv -f $PLATFORM $OUT/

done

echo Finished

#FQ1=$IN/${PLATFORM}_S1_L001_R1_001.fastq
#if [ -f "$FQ1" ]; then
#	FQ2=${FQ1/_R1_/_R2_}
#else
#	echo "$FQ1 does not exist"
#fi
#
## make a dummy read2 file
#cat $FQ1 | awk --posix '{ if (NR % 4 == 0) { sub(/.*/,"CCCCCCCC")} else if (NR % 2 == 0) { sub(/.*/,"NNNNNNNN")}; print}' > $FQ2
#
#
## IMPORTANT, because if dir already exists, longranger wont run again, and we might not notice
#rm -rf $PLATFORM       # where longranger will put it
#rm -rf $OUT/$PLATFORM  # where we will put it in the end
#
## Run longranger to demultplex single cells/gems
#echo "Running Longranger on $PLATFORM..."
#longranger basic --id=$PLATFORM --fastqs=$IN --PLATFORM=$PLATFORM > $OUT/$PLATFORM.log 2>&1
#
## Longranger places the files locally
#[[ -d $PLATFORM ]] && mv -f $PLATFORM $OUT/

