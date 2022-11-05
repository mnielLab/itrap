#!/bin/bash

#WRK_DIR=/home/tuba/herpov/tcr-pmhc-sc-project
#RAW_FIL=/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_IONTORRENT/raw/R_2019_10_29_17_32_41_user_GSS5PR-0198-13-4060_DTU_BSeq043_4060_DTU_BSeq043.fastq
#CUT_FIL=/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_IONTORRENT/processed/cut_adapter/IONTORRENT_S1_L001_R1_001.fastq
#LRO_FIL=/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_IONTORRENT/processed/longranger_out/IONTORRENT/outs/barcoded.fastq.gz
#LRC_FIL=/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_IONTORRENT/processed/longranger_clean/IONTORRENT.R1.gems.no_umi.fq
#KMA_FIL=/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_IONTORRENT/kma/kma-1t1/output/mapping.clean.gz
#
#CNT_FIL=/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_IONTORRENT/reports/read_counts_mhc_pipeline.txt

#WRK_DIR=/home/tuba/herpov/tcr-pmhc-sc-project
RAW_FIL=$1 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/raw/ILLUMINA-S1_L001_R1_001.fastq.gz
CUT_FIL=$2 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/cut_adapter/ILLUMINA_S1_L001_R1_001.fastq
LRO_FIL=$3 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/longranger_out/ILLUMINA/outs/barcoded.fastq.gz
LRC_FIL=$4 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/longranger_clean/ILLUMINA.R1.gems.no_umi.fq
KMA_FIL=$5 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/kma/kma-1t1/output/mapping.frag.gz

CNT_FIL=$6 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/reports/read_counts_mhc_pipeline.txt

echo -e "count\tfile" > $CNT_FIL

for FIL in $RAW_FIL $CUT_FIL $LRO_FIL $LRC_FIL $KMA_FIL
do
	var_name=$(dirname $FIL | rev | cut -d"/" -f1 | rev)
	# The "outs"-file also contains the faux read 2. Therefore we must divide by the double amount
	if [ $var_name = "outs" ]; then
		div=8
	else
		div=4
	fi

	if [ "${FIL#*.}" = "fastq.gz" ]; then
		count=$(echo $(zcat $FIL | wc -l )/$div|bc)
		echo -e $count"\t"$var_name >> $CNT_FIL

	elif [ "${FIL##*.}" = "fastq" ] || [ "${FIL##*.}" = "fq" ]; then
		count=$(echo $(cat $FIL | wc -l )/$div|bc)
		tcount=$()
		echo -e $count"\t"$var_name >> $CNT_FIL

	else
		count=$(zcat $FIL | wc -l )
		echo -e $count"\t"$var_name >> $CNT_FIL
		count=$(zcat $FIL | cut -f2 | awk '$1==1{c++} END{print c+0}')
		echo -e $count"\tcredible alignments" >> $CNT_FIL
	fi

	# Count reads longer than 100bp?
	#if [ $var_name = "outs" ]; then
	#	count=$(zcat $FIL | awk 'NR%4 == 2 {if (length > 100) print $0}' | wc -l)
	#	echo -e $count"\t100bp threshold" >> $CNT_FIL
	#fi
done



