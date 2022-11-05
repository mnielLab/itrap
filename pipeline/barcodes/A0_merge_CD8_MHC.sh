#EXP=exp2_MHC_ILLUMINA
#PARENT_DIR=/home/tuba/herpov/tcr-pmhc-sc-project/data/${EXP}/raw
MHC_IN=$1 #${PARENT_DIR}/MHC-2_S1_L001_R1_001.fastq.gz
CD8_IN=$2 #${PARENT_DIR}/CD8-2_S2_L001_R1_001.fastq.gz
FIL_OUT=$3 #${PARENT_DIR}/ILLUMINA-2_S1_L001_R1_001.fastq.gz

WRK_DIR=$(dirname MHC_IN)
MHC_OUT=${WRK_DIR}/MHC.fastq.gz
CD8_OUT=${WRK_DIR}/CD8.fastq.gz



zcat $MHC_IN | awk '{print (NR%4 == 1) ? $0 "-MHC" : $0}' | gzip -c > $MHC_OUT
zcat $CD8_IN | awk '{print (NR%4 == 1) ? $0 "-MHC" : $0}' | gzip -c > $CD8_OUT

cat $MHC_OUT $CD8_OUT > $FIL_OUT

rm $MHC_OUT $CD8_OUT
mkdir -p ${WRK_DIR}/raw
mv $MHC_IN $CD8_IN ${WRK_DIR}/raw/.

# OBS!
# EDIT script so that it concatenates all -fastq.gz files in a directory. The added string could be the basename of the file?