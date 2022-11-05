#!/usr/bin/env python
# coding: utf-8

from Bio import AlignIO
import pandas as pd
import numpy as np
import re

INPUT = snakemake.input[0] #'experiments/exp13/run1/mhc/mapping_split/align/A1072B317.aln'
OUTPUT = snakemake.output[0]

pattern = re.compile('N+')
matrix = list()

alignments = AlignIO.parse(INPUT, 'emboss')
for record in alignments:
    #print(record[1].id)
    #print(record[0].seq)
    #print(record[1].seq)
    match = pattern.search(str(record[0].seq))
    if match is None:
        umi_1 = ''
    else:
        idx1_s, idx1_e = match.span()
        umi_1 = str(record[1].seq)[idx1_s:idx1_e]    
        
    match = pattern.search(str(record[0].seq), idx1_e)
    if match is None:
        umi_2 = ''
    else:
        idx2_s, idx2_e = match.span()
        umi_2 = str(record[1].seq)[idx2_s:idx2_e]   
        
    match = re.search('-+' , umi_1 + umi_2)
    if match is None:
        umi = umi_1 + umi_2
    else:
        umi = ''
        
    matrix.append([record[0].id, record[1].id, umi])
    
pd.DataFrame(matrix).to_csv(OUTPUT, sep='\t', index=False, header=None)