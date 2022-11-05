#!/bin/bash
#set -euo pipefail
#IFS=$'\n\t'
# http://redsymbol.net/articles/unofficial-bash-strict-mode/#solution-positional-parameters

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

IN=$1
OUT=$2 #$PARENT_DIR/data/$EXP/processed/cut_adapter
PLATFORM=$3 #"ILLUMINA" #"IONTORRENT"

IN_DIR=$(dirname $IN)
OUT_DIR=$(dirname $OUT)

mkdir -p $OUT_DIR


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


#$PARENT_DIR/data/$EXP/raw/
for f in $IN_DIR/*fastq* ; do
	
# parentheses and the ampersand put this in a subshell, meaning that we can run all instances of longranger in parallel...!)

(	#---- Remove adapter sequences

	echo $f
	FN=$(basename $f)

	if [ "$PLATFORM" = "IONTORRENT" ]; then
		#FQ=$OUT_DIR/${PLATFORM}_$FILE_EXT #S1_L001_R1_001.fastq

		ADAPTER_5p="CACGACGCTCTTCCGATCT" #CACGACGCTCTTCCGATCT
		ADAPTER_3p="ATCACCGACTGCCCATAGAGAGG" #GGAGAGATACCCGTCAGCCACTA #CCTCTCTATGGGCAGTCGGTGAT

		#$CUTADAPT -a "^${ADAPTER_5p};optional...${ADAPTER_3p}" -o $OUT $f > $OUT_DIR/$PLATFORM.cutadapt.log # Added ;optional
		$CUTADAPT -a "^${ADAPTER_5p}...${ADAPTER_3p}" -o $OUT_DIR/$FN $f > $OUT_DIR/$FN.cutadapt.log # Added ;optional
	else
		# OBS! Instead of creating two files al the way down, merge the two together after annotating the read IDs with MHC or CD8.
		#FQ=$OUT_DIR/${PLATFORM}_$FILE_EXT #S1_L001_R1_001.fastq #-${BARCODE_TYPE}

		ADAPTER_3p="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"

		$CUTADAPT -a $ADAPTER_3p --error-rate 0.2 -o $OUT $f > $OUT_DIR/$PLATFORM.cutadapt.log
		#$TRIMGALORE --adapter $ADAPTER --quality 20 --length 60 --fastqc_args "--outdir $OUT/FastQC" --output_dir $OUT --path_to_cutadapt $CUTADAPT $f  > $OUT/$PLATFORM.trimgalore.log 2>&1
	fi
) &

done

wait

echo Finished
