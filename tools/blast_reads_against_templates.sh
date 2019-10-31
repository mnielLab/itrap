PLATFORM="ILLUMINA"
EXP="exp2_MHC_"$PLATFORM

PARENT_DIR=/home/tuba/herpov/tcr-pmhc-sc-project
#FASTA=$PARENT_DIR/data/$EXP/processed/longranger_clean/IONTORRENT.R1.gems.no_umi.no_adapters.fa
IN=$PARENT_DIR/data/$EXP/processed/longranger_clean
WRK_DIR=$PARENT_DIR/data/$EXP/blast #OBS! INVARIANT

#OVERREP_SEQS='TTTCTTATATGGGGCTCTCGGCCTTAGCGCCATTTTTTTGGAAACCTCTG|TTTCTTATATGGGCTCTCGGCCTTAGCGCCATTTTTTTGGAAACCTCTGC'


BLAST=/home/tuba-nobackup/shared/ncbi-blast/ncbi-blast-2.9.0+/bin

for sub_dir in $WRK_DIR/overrep_seq_templates $WRK_DIR/expected_templates #$WRK_DIR/rev_templates_start $WRK_DIR/reversed_templates
do
	echo $(basename $sub_dir)
	for fasta in $IN/*.fa #for fasta in $IN/ILLUMINA-CD8-2.R1.gems.no_umi.fa
	do
		echo $fasta
		BARCODE_TYPE=$(basename ${fasta%%.*})
		echo $BARCODE_TYPE

		$BLAST/makeblastdb -in $sub_dir/templates.fa -title "templates" -dbtype nucl
		# https://github.com/seqan/lambda/wiki/BLAST-Output-Formats
		$BLAST/blastn -db $sub_dir/templates.fa -query $fasta -outfmt 6 > $sub_dir/$BARCODE_TYPE.blast.tsv 2> $sub_dir/$BARCODE_TYPE.blast.err
		
		## Extract read IDs containing overrepresented sequences
		#cat $FASTA | paste - - | egrep $OVERREP_SEQS | cut -f1 | cut -f2 -d ">" | cut -f1 -d " " > $sub_dir/overrepresented_seq_IDs.lst
		#
		## Extract BLAST hits where query contains an overrepresented sequence
		#awk 'NR==FNR{arr[$1]++;next} arr[$1]' $sub_dir/overrepresented_seq_IDs.lst $sub_dir/blast.tsv > $sub_dir/overrep_in_blast_hits.tsv
	done
done

wait
