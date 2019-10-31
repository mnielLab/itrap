EXP=exp2_MHC_ILLUMINA
PARENT_DIR=/home/tuba/herpov/tcr-pmhc-sc-project/data/${EXP}/raw
MHC_IN=${PARENT_DIR}/MHC-2_S1_L001_R1_001.fastq.gz
MHC_OUT=${PARENT_DIR}/MHC.fastq.gz
CD8_IN=${PARENT_DIR}/CD8-2_S2_L001_R1_001.fastq.gz
CD8_OUT=${PARENT_DIR}/CD8.fastq.gz

FIL_OUT=${PARENT_DIR}/ILLUMINA-2_S1_L001_R1_001.fastq.gz

zcat $MHC_IN | awk '{print (NR%4 == 1) ? $0 "-MHC" : $0}' | gzip -c > $MHC_OUT
zcat $CD8_IN | awk '{print (NR%4 == 1) ? $0 "-MHC" : $0}' | gzip -c > $CD8_OUT

cat $MHC_OUT $CD8_OUT > $FIL_OUT

rm $MHC_OUT $CD8_OUT
mkdir -p raw
mv $MHC_IN $CD8_IN raw