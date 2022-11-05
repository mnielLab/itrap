#!/bin/bash

IN_DIR=$(echo $1 | rev | cut -d'/' -f3- | rev) #$(dirname $INP)
OUT_DIR=$(echo $2 | rev | cut -d'/' -f4- | rev)

echo $OUT_DIR

# Do I have to delete the folder or is it reminisance from testing stuf?
if [ -d "$OUT_DIR" ]; then
	echo "I exist"
	rm -rf $OUT_DIR
	mkdir $OUT_DIR
fi

cd $OUT_DIR

PRJ=$(find $IN_DIR/Reports/html/*  -maxdepth 0 -type d | sed 's/.*-//')
#REF=/home/tuba/marquard/tcr_seq/10x/apps/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0
#REF=/home/tuba-nobackup/shared/cellranger/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0
REF=/home/tuba-nobackup/shared/cellranger/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0

#CELLRANGER=/home/tuba-nobackup/shared/cellranger/cellranger-3.1.0/cellranger
#$CELLRANGER vdj --id=TCR_VDJ --fastqs=$IN_DIR/$PRJ/TCR/ --reference=$REF --sample=TCR --chain=TR

CELLRANGER=/home/tuba-nobackup/shared/cellranger/cellranger-6.1.1/cellranger
$CELLRANGER vdj --id=TCR_VDJ --fastqs=$IN_DIR/$PRJ/ --reference=$REF --sample=TCR --chain=TR
#$CELLRANGER vdj --id=TCR_VDJ_v4 --fastqs=$IN_DIR/$PRJ/TCR1/,$IN_DIR/$PRJ/TCR2/ --reference=$REF --sample=TCR1,TCR2 --chain=TR
#$CELLRANGER vdj --id=TCR_VDJ2 --fastqs=$IN_DIR/$PRJ/TCR2/ --reference=$REF --sample=TCR2 #--chain=TR

