



trim_galore --adapter TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG --quality 20 --length 0 --fastqc_args "--outdir tmp/FastQC" --output_dir tmp

/home/tuba/shared/bin/trim_galore --adapter TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG --quality 20 --length 0 --fastqc_args "--outdir tmp/FastQC" --output_dir tmp --path_to_cutadapt /home/people/marquard/.local/bin/cutadapt tmp/longranger_out/CD8-2/outs/barcoded.fastq





#!/bin/bash

# Get this script's parent dir, and work from there:
PARENT_DIR=`dirname $0`/../..


IN="$PARENT_DIR/tmp/longranger_clean"

# where to put tmp output of this script
OUT="$PARENT_DIR/tmp/barracoda_in"
mkdir -p $OUT

# where to put results output of this script
RESULTS="$PARENT_DIR/results"
mkdir -p $RESULTS



#-- OLIGO A B FASTA FILES
cp $PARENT_DIR/data/barcode_library/oligo_A.fa $RESULTS/a.fa
cp $PARENT_DIR/data/barcode_library/oligo_B.fa $RESULTS/b.fa

#-- ANNOTATIONS
cp $PARENT_DIR/data/barcode_library/CD8.annotations.xlsx $RESULTS/
cp $PARENT_DIR/data/barcode_library/MHC.annotations.xlsx $RESULTS/


IN_SUFFIX=".R1.gems.no_umi.fq"


for FILE in $IN/*${IN_SUFFIX}; do
# for SAMPLE in CD8 MHC; do
   
    SAMPLE=`basename $FILE ${IN_SUFFIX}`
    echo $FILE
    echo $SAMPLE


	ADAPTER=$(grep -v ">" $PARENT_DIR/data/barcode_library/$SAMPLE.seq-adapter.fa | $PARENT_DIR/scripts/utils/revcomp-seq.pl)

	# ---- Remove adapter sequences
	/home/people/marquard/.local/bin/cutadapt -a $ADAPTER -o $IN/$SAMPLE.R1.gems.no_umi.no_adapters.fq $FILE > $IN/$SAMPLE.cutadapt.log


	#---- Reverse complement all reads (and reverse quality line) -----------------------------

	$PARENT_DIR/scripts/utils/revcomp-fastq.pl "$IN/$SAMPLE.R1.gems.no_umi.no_adapters.fq" > "$OUT/$SAMPLE.R1.gems.no_umi.no_adapters.revcomp.fq"


	# add chromium barcode to start of read (no need to reverse complement it), so it becomes the sample ID for barracoda analysis
	# remember to add dummy chars to the quality line, too. The actual quality recorded for the chromium barcode nucleotides is lost, but that's no big deal.
	cat $OUT/$SAMPLE.R1.gems.no_umi.no_adapters.revcomp.fq | paste - - - - | sed -e 's/\([ACTG]\{16\}\)-1\t/\1-1\t\1/' | tr "\t" "\n" | awk --posix '{ if (NR % 4 == 0) { sub(/^/,"CCCCCCCCCCCCCCCC"); print} else {print}}' > $OUT/$SAMPLE.R1.gems.no_umi.no_adapters.revcomp.with_sample_id.fq


	#---- Get list of all chromium barcodes present and turn into sample fasta file, and samle key map for barracoda


	#-- SAMPLE ID FASTA
	# FASTA id is the same as the sequence AND the line number, to be compatible with barracoda
	cat $OUT/$SAMPLE.R1.gems.no_umi.no_adapters.revcomp.with_sample_id.fq | grep "BX:Z:" | sed -e 's/.*\([ACTG]\{16\}\).*/\1/' | sort | uniq | perl -ne 'chomp $_; print ">$_"."_$.\n"."$_\n";' > $RESULTS/$SAMPLE.sample_id.fasta 

	#-- SAMPLE KEY MAP
	cat  $RESULTS/$SAMPLE.sample_id.fasta  | paste - - | sed 's/>//' >  $RESULTS/$SAMPLE.sample-key-map.tsv
	# in this case 'sample' refers to chromium barcode / GEM / cell. I use 'sample' so it makes sense when we analyse these files with barracoda.




done
