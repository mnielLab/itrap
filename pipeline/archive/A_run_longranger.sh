#!/bin/bash

# Andrea Marquard
# 15.05.2019

# iontorrent reads start with part of the 10X-adapter Akey, which is:
#                    CCATCTCATCCCTGCGTGTCTCCGACTCAGCACGACGCTCTTCCGATC
# The reads begin like so:                         CACGACGCTCTTCCGATCT
# Then comes 26 bases (GEM barcode, and 16bp UMI)                     NNNNNNNNNNNNNNNNNNNNNNNNNN                       
# Then the TSO (# 4.5 M reads contain the 10x TSO)                                              TTTCTTATATGGG


# Get this script's parent dir, and work from there:
PARENT_DIR=`dirname $0`/../..

IN=$PARENT_DIR/tmp/longranger_in
OUT=$PARENT_DIR/tmp/longranger_out

mkdir -p $IN
mkdir -p $OUT


CUTADAPT="/home/people/marquard/.local/bin/cutadapt"
TRIMGALORE="/home/tuba/shared/bin/trim_galore"
    
ADAPTER="CACGACGCTCTTCCGATCT"




for f in $PARENT_DIR/data/fastq/*fastq ; do
	
# parentheses and the ampersand put this in a subshell, meaning that we can run all instances of longranger in parallel...!)

(	FILE=$(basename "${f}" .gz)
	# echo $FILE
	
	SAMPLE="IONTORRENT"
	# echo $SAMPLE

	# # gunzip fastq.gz files  # this data was not zipped..!
	FQ=$IN/IONTORRENT_S1_L001_R1_001.fastq 
	# gunzip -ck $f > $FQ


	# ---- Remove adapter sequences

    # $TRIMGALORE --adapter $ADAPTER --quality 20 --length 60 --fastqc_args "--outdir $IN/FastQC" --output_dir $IN --path_to_cutadapt $CUTADAPT $f  > $IN/$SAMPLE.trimgalore.log 2>&1
	# ..! trimgalore doesnt trim adaptors from 5' only 3'..!
	# ...So I use cutadapt instead:
	/home/people/marquard/.local/bin/cutadapt --front $ADAPTER --error-rate 0.2 -o $FQ $f > $IN/$SAMPLE.cutadapt.log


    # make a dummy read2 file
    FQ2=${FQ/_R1_/_R2_}
	cat $FQ | awk --posix '{ if (NR % 4 == 0) { sub(/.*/,"CCCCCCCC")} else if (NR % 2 == 0) { sub(/.*/,"NNNNNNNN")}; print}' > $FQ2


	# IMPORTANT, because if dir already exists, longranger wont run again, and we might not notice
	rm -rf $SAMPLE       # where longranger will put it
    rm -rf $OUT/$SAMPLE  # where we will put it in the end

	# Run longranger to demultplex single cells/gems
	echo "Running Longranger on $SAMPLE..."
	longranger basic --id=$SAMPLE --fastqs=$IN --sample=$SAMPLE > $OUT/$SAMPLE.log 2>&1

	[[ -d $SAMPLE ]] && mv -f $SAMPLE $OUT/

) &

done

wait


 #   1132 CACGACACTCTTCCGATCT
 #   1167 CACGGACGCTCTTCCGATCT
 #   1235 CACCGACGCTCTTCCGATCT
 #   1361 CACGACGCTTCCGATCT
 #   1430 CACGACGCTCTTCCG
 #   1505 CACGACGTTCTTCCGATCT
 #   1568 CACGACGCTTTTCCGATCT
 #   1575 CACGACGCTCTTTCGATCT
 #   1615 TACGACGCTCTTCCGATCT
 #   1838 CACGACGCTCTTCTGATCT
 #   1924 CACGAACGCTCTTCCGATCT
 #   1964 CACGATGCTCTTCCGATCT
 #   2403 CATGACGCTCTTCCGATCT
 #   2477 CAACGACGCTCTTCCGATCT
 #   2710 CACGACGGCTCTTCCGATCT
 #   3001 CACGACGCTCTTCCGA
 #   3223 CACACGCTCTTCCGATCT
 #   3435 CACGACCTCTTCCGATCT
 #   3600 CAGACGCTCTTCCGATCT
 #   3700 CACGACGCCTTCCGATCT
 #   4034 CACGACGTCTTCCGATCT
 #   4151 CCGACGCTCTTCCGATCT
 #   4591 CACGACGCTTTCCGATCT
 #   4612 CACGAGCTCTTCCGATCT
 #   4785 CACGCGCTCTTCCGATCT
 #   7789 CACGACGCTCTTCGATCT
 #   8571 ACGACGCTCTTCCGATCT
 #  10746 CACGACGCTCT CCGATCT
 #  12455 CACGACGCTCTTCCGAT
 # 120126 CACGACGCTCTTCCGATC
 #        CACGACGCTCTTCCGATCT
