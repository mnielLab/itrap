#!/bin/bash

KMA_FIL=$1 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/raw/ILLUMINA-S1_L001_R1_001.fastq.gz
SMP_FIL=$2 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/cut_adapter/ILLUMINA_S1_L001_R1_001.fastq
ALP_FIL=$3 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/longranger_out/ILLUMINA/outs/barcoded.fastq.gz
BET_FIL=$4 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/processed/longranger_clean/ILLUMINA.R1.gems.no_umi.fq
BCI_FIL=$5 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/kma/kma-1t1/output/mapping.frag.gz

UMI_FIL=$6 #/home/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/reports/read_counts_mhc_pipeline.txt

R_program=/home/tuba-nobackup/shared/R/R-3.6.1/bin/Rscript
parse_kma_results=/home/tuba/kamilla/10x-barcoding/scripts/parse-kma-results.R

$R_program ${parse_kma_results} $KMA_FIL $SMP_FIL $ALP_FIL $BET_FIL $BCI_FIL $UMI_FIL