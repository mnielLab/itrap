WRK_DIR="/home/tuba/herpov/tcr-pmhc-sc-project/data/exp3_MHC/kma_sandbox"

for kmer_size in 7 8 9 10 11 12 13
do
	for mismatch in -1 -2 -3 -4
	do
		for opening in 0 -1 -2 -3 -4 -5
		do
			for extension in 0 -1 -2 -3 -4 -5
			do
				ARGS="-1t1 -k $kmer_size -penalty $mismatch -gapopen $opening -gapextend $extension"
				KMA_DIR=$WRK_DIR"/kma"$(echo -e "${ARGS}" | tr -d '[:space:]')

				mapped_reads_count=$(zcat $KMA_DIR/output/mapping.frag.gz | wc -l | awk '{print $1}')

				echo -e "$ARGS\t$mapped_reads_count" >> $WRK_DIR/mapped_reads_count.lst
			done
		done
	done
done

sort -t$'\t' -k2 -n $WRK_DIR/mapped_reads_count.lst > $WRK_DIR/mapped_reads_count_sorted.lst