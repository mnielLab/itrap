#!/bin/bash

# Helle Povlsen
# 03.06.2019

PLATFORM="ILLUMINA"
EXP="exp2_MHC_"$PLATFORM # Variable
PROJECT=""
#BARCODE_TYPE=${PLATFORM}-MHC-2

PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
IN_DIR=$PARENT_DIR/data/$EXP/processed/longranger_clean
OUT_DIR=$PARENT_DIR/data/$EXP/cutadapt_naming #/$PROJECT

mkdir -p $OUT_DIR



CUTADAPT="/opt/anaconda3/bin/cutadapt"

# sequences:
TSO="TTTCTTATATGGG"
PRIMER_B="CGAGTACCATGGGCGTAAC"
ANNEAL="CTACGCCTTTTGGGGAAGGTCACAC"
PRIMER_MHC="GACGCTGGCTGGAACTTC"
PRIMER_CD8="GATGAAACCTGGGATCTTTTAC"

# OBS! Introduce adapter trimming here? Thus, quality trimmed reads will have adapters removed. According to manual, adapter trimming must be done after quality trimming.
for BARCODE_TYPE in ${PLATFORM}-MHC-2 ${PLATFORM}-CD8-2
do

	SUFFIX=$BARCODE_TYPE.R1.gems.no_umi #.$PROJECT
	# ---- Name reads
	echo "TSO"
	$CUTADAPT -g TSO=$TSO -e 0.1 -y ' {name}' --action=lowercase -o $OUT_DIR/${SUFFIX}.TSO.fq $IN_DIR/${SUFFIX}.fq > $OUT_DIR/cutadapt_naming.log
	echo "PRIMER_B"
	$CUTADAPT -g PRIMER_B=$PRIMER_B -e 0.1 -y ' {name}' --action=lowercase -o $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.fq $OUT_DIR/${SUFFIX}.TSO.fq >> $OUT_DIR/cutadapt_naming.log
	echo "ANNEAL"
	$CUTADAPT -g ANNEAL=$ANNEAL -e 0.1 -y ' {name}' --action=lowercase -o $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.fq $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.fq >> $OUT_DIR/cutadapt_naming.log
	echo "PRIMER_CD8"
	$CUTADAPT -a PRIMER_CD8=$PRIMER_CD8 -e 0.1 -y ' {name}' --action=lowercase -o $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8.fq $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.fq >> $OUT_DIR/cutadapt_naming.log
	echo "PRIMER_MHC"
	$CUTADAPT -a PRIMER_MHC=$PRIMER_MHC -e 0.1 -y ' {name}' --action=lowercase -o $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.fq $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8.fq >> $OUT_DIR/cutadapt_naming.log

	sed -n '1~4s/^@/>/p;2~4p' $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.fq > $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.fa

	# Making headers
	mkdir -p $OUT_DIR/headers
	grep "HBBAF:" $OUT_DIR/${SUFFIX}.TSO.PRIMER_B.ANNEAL.PRIMER_CD8_MHC.fq | cut -d"@" -f2 > $OUT_DIR/headers/${BARCODE_TYPE}.annotated_headers.lst
	cut -d" "  -f3- $OUT_DIR/headers/annotated_headers.lst | sort | uniq -c | sed 's/^ *//g' > $OUT_DIR/headers/${BARCODE_TYPE}.annotated_headers.uniq_count.lst

done
