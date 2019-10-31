#!/bin/bash

# Andrea Marquard
# 15.05.2019

PROJECT=IONTORRENT

# Get this script's parent dir, and work from there:
PARENT_DIR=`dirname $0`/../..


IN="$PARENT_DIR/tmp/longranger_clean"

# where to put tmp output of this script
OUT="$PARENT_DIR/tmp/barracoda_in"
mkdir -p $OUT

# where to put results output of this script
RESULTS="$PARENT_DIR/results"
mkdir -p $RESULTS


IN_SUFFIX=".R1.gems.no_umi"

mkdir -p $IN/FastQC

CUTADAPT="/home/people/marquard/.local/bin/cutadapt"

FILE=$IN/$PROJECT${IN_SUFFIX}.fq



# MHC_P1_ID1 CCTCTCTATGGGCAGTCGGTGATCTCTCCAAGAAGTTCCAGCCAGCGTC
# CD8_P1_ID5 CCTCTCTATGGGCAGTCGGTGATAGGCTGAAGTAAAAGATCCCAGGTTTCATC
#            |||||||||||||||||||||||
#            RevComp of this is the adapter to remove
ADAPTER=ATCACCGACTGCCCATAGAGAGG


# ---- Remove adapter sequences

$CUTADAPT -a $ADAPTER -o $IN/$PROJECT${IN_SUFFIX}.no_adapters.fq $FILE > $IN/$PROJECT.cutadapt.log



#---- Reverse complement all reads (and reverse quality line) -----------------------------

$PARENT_DIR/scripts/utils/revcomp-fastq.pl "$IN/$PROJECT${IN_SUFFIX}.no_adapters.fq" > "$OUT/$PROJECT${IN_SUFFIX}.no_adapters.revcomp.fq"


