#!/bin/bash

# Helle Povlsen
# 02.06.19

# iontorrent reads start with part of the 10X-adapter Akey, which is:
#                    CCATCTCATCCCTGCGTGTCTCCGACTCAGCACGACGCTCTTCCGATC
# The reads begin like so with 5' adapter:         CACGACGCTCTTCCGATCT
# Then comes 26 bases (GEM barcode, and 16bp UMI)                     NNNNNNNNNNNNNNNNNNNNNNNNNN                       
# Then the TSO (# 4.5 M reads contain the 10x TSO)                                              TTTCTTATATGGG

# Reads end with 3' adapter
# MHC_P1_ID1 CCTCTCTATGGGCAGTCGGTGATCTCTCCAAGAAGTTCCAGCCAGCGTC
# CD8_P1_ID5 CCTCTCTATGGGCAGTCGGTGATAGGCTGAAGTAAAAGATCCCAGGTTTCATC
#            |||||||||||||||||||||||
#            RevComp of this is the 3' adapter to remove

PLATFORM="ILLUMINA" #"IONTORRENT"
EXP="exp2_MHC_"$PLATFORM # Variable

PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
OUT=$PARENT_DIR/data/$EXP/processed/cut_adapter

mkdir -p $OUT


# This script uses cutadapt to trim adaptor sequences, since trimaglore does not trim adaptors from 5', but only from 3'
CUTADAPT="/opt/anaconda3/bin/cutadapt"
TRIMGALORE="/home/tuba/shared/bin/trim_galore"

# https://cutadapt.readthedocs.io/en/stable/guide.html
# Cutadapt options
# --quality-cutoff 20 (before adapter removal)
#Filtering of processed reads:
#  Filters are applied after above read modifications.
#
#  -m LEN[:LEN2], --minimum-length LEN[:LEN2]
#                        Discard reads shorter than LEN. Default: 0
#  -M LEN[:LEN2], --maximum-length LEN[:LEN2]
#                        Discard reads longer than LEN. Default: no limit
#
# Output
#  --too-short-output FILE
#                        Write reads that are too short (according to length
#                        specified by -m) to FILE. Default: discard reads
#  --too-long-output FILE
#                        Write reads that are too long (according to length
#                        specified by -M) to FILE. Default: discard reads

    


#FQ1=$OUT/${PLATFORM}_no_5p_adapter.fastq 
#FQ2=$OUT/${PLATFORM}_S1_L001_R1_001.fastq
#FQ=$OUT/${PLATFORM}_S1_L001_R1_001.fastq

for f in $PARENT_DIR/data/$EXP/raw/*fastq* ; do
	
# parentheses and the ampersand put this in a subshell, meaning that we can run all instances of longranger in parallel...!)

(	#$CUTADAPT --front $ADAPTER_5p --error-rate 0.2 -o $FQ1 $f > $OUT/$PLATFORM.cutadapt_5p.log
	#$CUTADAPT -a $ADAPTER_3p -o $FQ2 $FQ1 > $OUT/$PLATFORM.cutadapt_3p.log

	#CUTADAPT -a ${ADAPTER_5p}...${ADAPTER_3p} -o FQ2 $f > $OUT/$PLATFORM.cutadapt.log

	#---- Remove adapter sequences

	if [ "$PLATFORM" = "IONTORRENT" ]; then
		FQ=$OUT/${PLATFORM}_S1_L001_R1_001.fastq

		ADAPTER_5p="CACGACGCTCTTCCGATCT"
		ADAPTER_3p="ATCACCGACTGCCCATAGAGAGG"

		CUTADAPT -a ${ADAPTER_5p}...${ADAPTER_3p} -o FQ $f > $OUT/$PLATFORM.cutadapt.log
	else
		# OBS! Instead of creating two files al the way down, merge the two together after annotating the read IDs with MHC or CD8.
		#FILE=$(basename "${f}")
		#BARCODE_TYPE=${FILE/_*/}
		#echo $BARCODE_TYPE
		FQ=$OUT/${PLATFORM}_S1_L001_R1_001.fastq #-${BARCODE_TYPE}

		ADAPTER_3p="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"

		$CUTADAPT -a $ADAPTER_3p --error-rate 0.2 -o $FQ $f > $OUT/$PLATFORM.cutadapt.log
		#$TRIMGALORE --adapter $ADAPTER --quality 20 --length 60 --fastqc_args "--outdir $OUT/FastQC" --output_dir $OUT --path_to_cutadapt $CUTADAPT $f  > $OUT/$PLATFORM.trimgalore.log 2>&1
	fi
) &

done

wait

echo Finished
