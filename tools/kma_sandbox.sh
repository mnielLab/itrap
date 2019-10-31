for kmer_size in 7 8 9 10 11 12 13
do
	for mismatch in -1 -2 -3 -4
	do
		for opening in 0 -1 -2 -3 -4 -5
		do
			for extension in 0 -1 -2 -3 -4 -5
			do
				ARGS="-1t1 -k $kmer_size -penalty $mismatch -gapopen $opening -gapextend $extension"
				KMA_DIR="/home/tuba/herpov/tcr-pmhc-sc-project/data/exp3_MHC/kma_sandbox/kma"$(echo -e "${ARGS}" | tr -d '[:space:]')
				BARCODE_TEMPLATES="/home/tuba/herpov/tcr-pmhc-sc-project/data/exp3_MHC/barcode_library/barcode_templates.fa"
				READS="/home/tuba/herpov/tcr-pmhc-sc-project/data/exp3_MHC/processed/tmp/sandbox_reads.fq"

				KMA=/home/tuba-nobackup/shared/src/kma/kma

				echo $ARGS
				mkdir -p $KMA_DIR/output

				$KMA index -i $BARCODE_TEMPLATES -o $KMA_DIR/templates
				$KMA -i $READS -o $KMA_DIR/output/mapping -t_db $KMA_DIR/templates $ARGS

			done
		done
	done
done

# ended the program at kma-1t1-k7-penalty-4-gapopen-4-gapextend-3 (resume from here) --> reinsert kmer size 7
