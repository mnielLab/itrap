#!/bin/bash

# Andrea Marquard
# 15.05.2019

# Run fastqc to see how much adapter contamination there is 
# (I've seen some small RNA RPI1 adapter)


# Get this script's parent dir, and work from there:
PARENT_DIR=`dirname $0`/../..

# IN=$PARENT_DIR/tmp/longranger_in
IN=$PARENT_DIR/tmp/longranger_clean

# where to put output from this script
OUT="$PARENT_DIR/results/fastqc/"
mkdir -p $OUT



for i in $IN/*R1*.fq; do
	
(	fastqc --extract -o $OUT $i )&

done

wait
