#!/bin/bash

#WRK_DIR=/home/tuba/herpov/tcr-pmhc-sc-project
RAW_FIL=$1 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/raw/ILLUMINA-S1_L001_R1_001.fastq.gz
CUT_FIL=$2 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/cut_adapter/ILLUMINA_S1_L001_R1_001.fastq
#LRO_FIL=$3 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/longranger_out/ILLUMINA/outs/barcoded.fastq.gz
LRC_FIL=$3 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/longranger_clean/ILLUMINA.R1.gems.no_umi.fq
KMA_FIL=$4 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/kma/kma-1t1/output/mapping.frag.gz

OUT_DIR=$(echo "${@: -1}") #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/reports/read_counts_mhc_pipeline.txt


for FIL in $RAW_FIL $CUT_FIL $LRC_FIL $KMA_FIL
do
	var_name=$(dirname $FIL | rev | cut -d"/" -f1 | rev)
	# The "outs"-file also contains the faux read 2. Therefore we must divide by the double amount
	if [ "${FIL#*.}" = "fastq.gz" ]; then
		zcat $FIL | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $OUT_DIR/$var_name.txt

	elif [ "${FIL##*.}" = "fastq" ] || [ "${FIL##*.}" = "fq" ]; then
		cat $FIL | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $OUT_DIR/$var_name.txt

	else
		zcat $FIL | cut -f1 | awk '{print length($1)}' | sort -n | uniq -c > $OUT_DIR/$var_name.txt
	fi
done



