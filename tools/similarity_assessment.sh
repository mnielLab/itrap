#!/bin/bash

# Run seq2score
seq2score=/home/tuba-nobackup/shared/seq2score_db_kernel 
BLF=/home/tuba-nobackup/shared/BLOSUM50
QIJ=/home/tuba-nobackup/shared/blosum62.qij
F1=lol.txt
F2=lol.txt
F_OUT=lol.results.txt

$seq2score -blf $BLF -blqij $QIJ -pa $F1 $F2 > $F_OUT
