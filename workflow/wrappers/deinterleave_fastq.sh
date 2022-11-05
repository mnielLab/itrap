interleaved_fq=$1
read1_fq=$2
read2_fq=$3

zcat $interleaved_fq | paste - - - - - - - - | grep "BX:Z:" \
    | tee >(cut -f 1-4 | tr "\t" "\n" | gzip > $read1_fq) \
    |       cut -f 5-8 | tr "\t" "\n" | gzip > $read2_fq
