#!/bin/bash

# Helle Povlsen
# 03.06.2019

EXP="exp2_MHC" # Variable

PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
DIR=$PARENT_DIR/data/$EXP/processed/longranger_clean

PROJECT=IONTORRENT
IN_SUFFIX=".R1.gems.no_umi"
FILE=$DIR/$PROJECT${IN_SUFFIX}.fq

CUTADAPT="/opt/anaconda3/bin/cutadapt"

# MHC_P1_ID1 CCTCTCTATGGGCAGTCGGTGATCTCTCCAAGAAGTTCCAGCCAGCGTC
# CD8_P1_ID5 CCTCTCTATGGGCAGTCGGTGATAGGCTGAAGTAAAAGATCCCAGGTTTCATC
#            |||||||||||||||||||||||
#            RevComp of this is the adapter to remove
ADAPTER=ATCACCGACTGCCCATAGAGAGG


# ---- Remove adapter sequences

$CUTADAPT -a $ADAPTER -o $DIR/$PROJECT${IN_SUFFIX}.no_adapters.fq $FILE > $DIR/$PROJECT.cutadapt.log

sed -n '1~4s/^@/>/p;2~4p' $DIR/$PROJECT${IN_SUFFIX}.no_adapters.fq > $DIR/$PROJECT${IN_SUFFIX}.no_adapters.fa



#/opt/anaconda3/bin/cutadapt -q 20 -m 113 -M 155 --too-short-output IONTORRENT.R1.gems.no_umi.quality_trimmed.too_short.fq --too-long-output IONTORRENT.R1.gems.no_umi.quality_trimmed.too_long.fq -o IONTORRENT.R1.gems.no_umi.quality_trimmed.fq IONTORRENT.R1.gems.no_umi.fq
# Do not remove sequences too long? These sequences still contain adapters! Therefore they may be a bit longer than what I have anticipated here..