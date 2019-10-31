#!/bin/bash

# Helle Povlsen
# 08.07.2019

# Make stats for barcode building blocks

EXP="exp3_MHC" # Variable

PARENT_DIR="/home/tuba/herpov/tcr-pmhc-sc-project"
IN_MHC=$PARENT_DIR/data/$EXP/
OUT=$PARENT_DIR/data/$EXP/barcode_buildingblock_stats

PRIMER_10X="CACGACGCTCTTCCGATCT"
TSO=TTTCTTATATGGG
PRIMER_B=CGAGTACCATGGGCGTAAC
OLIGO_B='CTTGGCAATCCATGCTCCCATTTGG|GAACCATGTGGCTTGCTTCAACTTC|GCCTGTAGTCCCACGCGATCTAACA|GTTTACGTGAGTTGGAGGACTTTAT|GTATCAAGAGACGCTCTACCGGGAT|TGTGTCGCTGAAAGAACTCATGTCG|CAGTAAGTCGCTCCGGGTGGACTTG|CAATTACAGGGCAGGCAGGGTTCTC|GTTTCTTGGCCTCCAGCTAGTACAG|GTACTAATACTGAACGATCAGCCGC|GCTATCCTTGAGGCGGTGTGTACGG|CGATTAGAAGAGTCTCCGTTCACGT|TTTGCCGTTTCTGAGCCACTAATGA|GGGAACTATCCCTATTTTCAAACTC|GGCTGACCGGGAGGCTTGAATTGTA'
ANNEAL=CTACGCCTTTTGGGGAAGGTCACAC
OLIGO_A='TTGTCATAAGGAGATAGCTACTACG|CGGTAGTTACTTGCACTTTGCGGTC|GCTTGCAGGCAGATAATAACAAGCG|GTATCATGCAGCGGAATCTGTGCGC|ACAACCAGGTATACCCTCTTGGACC|TGGCTCCCCAGTCATTAAGAAACGT|CGGGAAACACCATCCAATTAGCCGC|CAATGCCCAGTTATGGAGTAACCTC|GGTTATCAAGGTCCTATACCCGCAG'
PRIMER_A_MHC=GACGCTGGCTGGAACTTC
PRIMER_A_CD8=GATGAAACCTGGGATCTTTTAC
SAMPLE_ID='TTGGAGAG|GCCGTTAA|TTGCTTGA|GCACAGGA|TTCAGCCT'
ADAPTER_3P=ATCACCGACTGCCCATAGAGAGG

OVERREP_SEQ1=TTTCTTATATGGGGCTCTCGGCCTTAGCGCCATTTTTTTGGAAACCTCTG
OVERREP_SEQ2=TTTCTTATATGGGCTCTCGGCCTTAGCGCCATTTTTTTGGAAACCTCTGC
OVERREP_SEQ3=TTTCTTATATGGGGCCTTTTTGCTGTAGGCCCGGGTGGTTGCTGCCGAAA


for file in `find $IN_MHC -type f \( -iname "*.fastq" -o -iname "*.fq" \)`; do

(	DIR="$OUT/$(basename ${file} .${file##*.})"
	mkdir -p $DIR

	echo $file
	STAT_FILE=${DIR}/stat_file.tab

	# Get read lengths
	grep "^"$PRIMER_10X $file | awk '{ print length($0) }' > ${DIR}/read_lengths_primer10x.lst
	grep "^"$TSO $file | awk '{ print length($0) }' > ${DIR}/read_lengths_TSO.lst
	grep $PRIMER_B $file | awk '{ print length($0) }' > ${DIR}/read_lengths_primerB.lst
	grep "^"$TSO$PRIMER_B $file | awk '{ print length($0) }' > ${DIR}/read_lengths_TSO+primerB.lst
	grep $ANNEAL $file | awk '{ print length($0) }' > ${DIR}/read_lengths_anneal.lst
	grep $PRIMER_A_MHC $file | awk '{ print length($0) }' > ${DIR}/read_lengths_primerMHC.lst
	grep $PRIMER_A_CD8 $file | awk '{ print length($0) }' > ${DIR}/read_lengths_primerCD8.lst
	grep $ADAPTER_3P $file | awk '{ print length($0) }' > ${DIR}/read_lengths_adapter_3p.lst
	grep $OVERREP_SEQ1 $file | awk '{ print length($0) }' > ${DIR}/read_lengths_overrep_seq1.lst
	grep $OVERREP_SEQ2 $file | awk '{ print length($0) }' > ${DIR}/read_lengths_overrep_seq2.lst
	grep $OVERREP_SEQ3 $file | awk '{ print length($0) }' > ${DIR}/read_lengths_overrep_seq3.lst

	egrep $OLIGO_B $file | awk '{ print length($0) }' > ${DIR}/read_lengths_oligoB.lst
	egrep $OLIGO_A $file | awk '{ print length($0) }' > ${DIR}/read_lengths_oligoA.lst
	egrep $SAMPLE_ID $file | awk '{ print length($0) }' > ${DIR}/read_lengths_sampleid.lst
	egrep $OLIGO_B $file | egrep $OLIGO_A | awk '{ print length($0) }' > ${DIR}/read_lengths_oligoAB.lst

	# Get stats in the form of line counts alias the number of reads
	n_total=$(($(wc -l $file | awk '{print $1}') / 4))
	n_10x_primer=$(wc -l ${DIR}/read_lengths_primer10x.lst | awk '{print $1}')
	n_tso=$(wc -l ${DIR}/read_lengths_TSO.lst | awk '{print $1}')
	n_primer_b=$(wc -l ${DIR}/read_lengths_primerB.lst | awk '{print $1}')
	n_tso_primer_b==$(wc -l ${DIR}/read_lengths_TSO+primerB.lst | awk '{print $1}')
	n_oligo_b=$(wc -l ${DIR}/read_lengths_oligoB.lst | awk '{print $1}')
	n_anneal=$(wc -l ${DIR}/read_lengths_anneal.lst | awk '{print $1}')
	n_oligo_a=$(wc -l ${DIR}/read_lengths_oligoA.lst | awk '{print $1}')
	n_primer_mhc=$(wc -l ${DIR}/read_lengths_primerMHC.lst | awk '{print $1}')
	n_primer_cd8=$(wc -l ${DIR}/read_lengths_primerCD8.lst | awk '{print $1}')
	n_sample_id=$(wc -l ${DIR}/read_lengths_sampleid.lst | awk '{print $1}')
	n_adapter_3p=$(wc -l ${DIR}/read_lengths_adapter_3p.lst | awk '{print $1}')
	n_oligo_ab=$(wc -l ${DIR}/read_lengths_oligoAB.lst | awk '{print $1}')
	n_overrep_seq1=$(wc -l ${DIR}/read_lengths_overrep_seq1.lst | awk '{print $1}')
	n_overrep_seq2=$(wc -l ${DIR}/read_lengths_overrep_seq2.lst | awk '{print $1}')
	n_overrep_seq3=$(wc -l ${DIR}/read_lengths_overrep_seq3.lst | awk '{print $1}')

	# Write building block stats
	echo -e "filter\tseq\t$(basename ${file} .${file##*.})" > $STAT_FILE
	echo -e "total\t-\t$n_total" >> $STAT_FILE
	echo -e "PRIMER_10X\t$PRIMER_10X\t$n_10x_primer" >> $STAT_FILE
	echo -e "TSO\t$TSO\t$n_tso" >> $STAT_FILE
	echo -e "PRIMER_B\t$PRIMER_B\t$n_primer_b" >> $STAT_FILE
	echo -e "TSO+PRIMER_B\t$TSO$PRIMER_B\t$n_tso_primer_b" >> $STAT_FILE
	echo -e "OLIGO_B\t-\t$n_oligo_b" >> $STAT_FILE
	echo -e "ANNEAL\t$ANNEAL\t$n_anneal" >> $STAT_FILE
	echo -e "OLIGO_A\t-\t$n_oligo_a" >> $STAT_FILE
	echo -e "PRIMER_A_MHC\t$PRIMER_A_MHC\t$n_primer_mhc" >> $STAT_FILE
	echo -e "PRIMER_A_CD8\t$PRIMER_A_CD8\t$n_primer_cd8" >> $STAT_FILE
	echo -e "SAMPLE_ID\t-\t$n_sample_id" >> $STAT_FILE
	echo -e "OLIGO_AB\t-\t$n_oligo_ab" >> $STAT_FILE
	echo -e "ADAPTER_3P\t$ADAPTER_3P\t$n_adapter_3p" >> $STAT_FILE
	echo -e "OVERREP_SEQ1\t$OVERREP_SEQ1\t$n_overrep_seq1" >> $STAT_FILE
	echo -e "OVERREP_SEQ2\t$OVERREP_SEQ2\t$n_overrep_seq2" >> $STAT_FILE
	echo -e "OVERREP_SEQ3\t$OVERREP_SEQ3\t$n_overrep_seq3" >> $STAT_FILE

)&

done

IN_TCR=${IN_MHC/_MHC/_TCR}

for file in `find $IN_TCR -type f \( -iname "*.fastq.gz" -o -iname "*.fq" \)`; do

(	DIR="$OUT/$(basename ${file} .${file##*.})"
	mkdir -p $DIR

	echo $file
	STAT_FILE=${DIR}/stat_file.tab

	# Get read lengths
	zcat $file | grep "^"$PRIMER_10X | awk '{ print length($0) }' > ${DIR}/read_lengths_primer10x.lst
	zcat $file | grep "^"$TSO | awk '{ print length($0) }' > ${DIR}/read_lengths_TSO.lst
	zcat $file | grep $PRIMER_B | awk '{ print length($0) }' > ${DIR}/read_lengths_primerB.lst
	zcat $file | grep "^"$TSO$PRIMER_B | awk '{ print length($0) }' > ${DIR}/read_lengths_TSO+primerB.lst
	zcat $file | grep $ANNEAL | awk '{ print length($0) }' > ${DIR}/read_lengths_anneal.lst
	zcat $file | grep $PRIMER_A_MHC | awk '{ print length($0) }' > ${DIR}/read_lengths_primerMHC.lst
	zcat $file | grep $PRIMER_A_CD8 | awk '{ print length($0) }' > ${DIR}/read_lengths_primerCD8.lst
	zcat $file | grep $ADAPTER_3P | awk '{ print length($0) }' > ${DIR}/read_lengths_adapter_3p.lst
	zcat $file | grep $OVERREP_SEQ1 | awk '{ print length($0) }' > ${DIR}/read_lengths_overrep_seq1.lst
	zcat $file | grep $OVERREP_SEQ2 | awk '{ print length($0) }' > ${DIR}/read_lengths_overrep_seq2.lst
	zcat $file | grep $OVERREP_SEQ3 | awk '{ print length($0) }' > ${DIR}/read_lengths_overrep_seq3.lst

	zcat $file | egrep $OLIGO_B | awk '{ print length($0) }' > ${DIR}/read_lengths_oligoB.lst
	zcat $file | egrep $OLIGO_A | awk '{ print length($0) }' > ${DIR}/read_lengths_oligoA.lst
	zcat $file | egrep $SAMPLE_ID | awk '{ print length($0) }' > ${DIR}/read_lengths_sampleid.lst
	zcat $file | egrep $OLIGO_B | egrep $OLIGO_A | awk '{ print length($0) }' > ${DIR}/read_lengths_oligoAB.lst

	# Get stats in the form of line counts alias the number of reads
	n_total=$(($(zcat $file | wc -l | awk '{print $1}') / 4))
	n_10x_primer=$(wc -l ${DIR}/read_lengths_primer10x.lst | awk '{print $1}')
	n_tso=$(wc -l ${DIR}/read_lengths_TSO.lst | awk '{print $1}')
	n_primer_b=$(wc -l ${DIR}/read_lengths_primerB.lst | awk '{print $1}')
	n_tso_primer_b=$(wc -l ${DIR}/read_lengths_TSO+primerB.lst | awk '{print $1}')
	n_oligo_b=$(wc -l ${DIR}/read_lengths_oligoB.lst | awk '{print $1}')
	n_anneal=$(wc -l ${DIR}/read_lengths_anneal.lst | awk '{print $1}')
	n_oligo_a=$(wc -l ${DIR}/read_lengths_oligoA.lst | awk '{print $1}')
	n_primer_mhc=$(wc -l ${DIR}/read_lengths_primerMHC.lst | awk '{print $1}')
	n_primer_cd8=$(wc -l ${DIR}/read_lengths_primerCD8.lst | awk '{print $1}')
	n_sample_id=$(wc -l ${DIR}/read_lengths_sampleid.lst | awk '{print $1}')
	n_adapter_3p=$(wc -l ${DIR}/read_lengths_adapter_3p.lst | awk '{print $1}')
	n_oligo_ab=$(wc -l ${DIR}/read_lengths_oligoAB.lst | awk '{print $1}')
	n_overrep_seq1=$(wc -l ${DIR}/read_lengths_overrep_seq1.lst | awk '{print $1}')
	n_overrep_seq2=$(wc -l ${DIR}/read_lengths_overrep_seq2.lst | awk '{print $1}')
	n_overrep_seq3=$(wc -l ${DIR}/read_lengths_overrep_seq3.lst | awk '{print $1}')

	# Write building block stats
	echo -e "filter\tseq\t$(basename ${file} .${file##*.})" > $STAT_FILE
	echo -e "total\t-\t$n_total" >> $STAT_FILE
	echo -e "PRIMER_10X\t$PRIMER_10X\t$n_10x_primer" >> $STAT_FILE
	echo -e "TSO\t$TSO\t$n_tso" >> $STAT_FILE
	echo -e "PRIMER_B\t$PRIMER_B\t$n_primer_b" >> $STAT_FILE
	echo -e "TSO+PRIMER_B\t$TSO$PRIMER_B\t$n_tso_primer_b" >> $STAT_FILE
	echo -e "OLIGO_B\t-\t$n_oligo_b" >> $STAT_FILE
	echo -e "ANNEAL\t$ANNEAL\t$n_anneal" >> $STAT_FILE
	echo -e "OLIGO_A\t-\t$n_oligo_a" >> $STAT_FILE
	echo -e "PRIMER_A_MHC\t$PRIMER_A_MHC\t$n_primer_mhc" >> $STAT_FILE
	echo -e "PRIMER_A_CD8\t$PRIMER_A_CD8\t$n_primer_cd8" >> $STAT_FILE
	echo -e "SAMPLE_ID\t-\t$n_sample_id" >> $STAT_FILE
	echo -e "OLIGO_AB\t-\t$n_oligo_ab" >> $STAT_FILE
	echo -e "ADAPTER_3P\t$ADAPTER_3P\t$n_adapter_3p" >> $STAT_FILE
	echo -e "OVERREP_SEQ1\t$OVERREP_SEQ1\t$n_overrep_seq1" >> $STAT_FILE
	echo -e "OVERREP_SEQ2\t$OVERREP_SEQ2\t$n_overrep_seq2" >> $STAT_FILE
	echo -e "OVERREP_SEQ3\t$OVERREP_SEQ3\t$n_overrep_seq3" >> $STAT_FILE

)&

done

wait

