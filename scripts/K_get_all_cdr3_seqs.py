#!/usr/bin/env python

import numpy as np
import pandas as pd
import re

def cdr3_lst_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace("'","").split(" ")

converters={'cdr3_lst_TRA': cdr3_lst_converter, 'cdr3_lst_TRB': cdr3_lst_converter}



def evaluate_cdr3(seqs):
    regex = re.compile('^[GPAVLIMCFYWHKRQNEDST]+$')
    
    cdr3s = list()
    for seq in seqs:
        match = regex.search(seq)
        
        if match is not None:
            cdr3s.append(seq)
    return np.unique(cdr3s)
        

# Input
IN_FILE = snakemake.params[0]
CHAIN = snakemake.params.chain
OUTPUT = snakemake.output[0]

# Load data
df = pd.read_csv(IN_FILE, usecols=['chain','cdr3'])

# OBS! The CDR3 sequences CANNOT contain characters other than the 20 AA letters!
if CHAIN.upper() in 'TRA':
    cdr3s = evaluate_cdr3(df[df.chain == 'TRA'].cdr3.dropna())
elif CHAIN.upper() in 'TRB':
    cdr3s = evaluate_cdr3(df[df.chain == 'TRB'].cdr3.dropna())
else:
    print('chain doesnt match anything')

np.savetxt(OUTPUT, cdr3s, fmt='%s')
#np.savetxt('all_cdr3_a', TRAs, fmt='%s') # permanently save, tmp_A2
#np.savetxt('all_cdr3_b', TRBs, fmt='%s') # permanently save

#####################
# END OF SCRIPT ONE #
#####################
        
# # Write chunck files for parallel processing
# for chain, lst in [('A', TRAs), ('B', TRBs)]:
#     # Get chunks
#     n = 50
#     for i in range(0, len(lst), n):
#         np.savetxt(f'{chain}_chunk{i}', lst[i:i + n], fmt='%s') # permanently save
# 
# 
# 
# ######################
# #     SCRIPT TWO     #
# ######################
# 
# # Run in parallel: Chains and chunks
# 
# import numpy as np
# import pandas as pd
# import subprocess
# import os
# 
# def run_cmd(cmd, input_string=''):
#     """Run the cmd with input_string as stdin and return output."""
#     p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
#                                              stderr=subprocess.PIPE, universal_newlines=True, close_fds=True)
#     result = pd.read_csv(p.stdout, sep=" ", names=['seq1', 'seq2', 'similarity'], usecols=[1,2,3], comment='#')
#     p.wait()
#     return result
# 
# # Input
# seq2score = snakemake.input.s2s #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel" #
# BLF = snakemake.input.blf #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50" #
# QIJ = snakemake.input.qij #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij" #
# 
# CHAIN = snakemake.params.chain
# CHUNK_I = snakemake.params.chunk_i
# CHUNK = snakemake.input.chu
# ALL_CDR3 = snakemake.input.cdr
# 
# OUTPUT = snakemake.output.chunk
# 
# chunk = np.loadtxt(CHUNK, fmt='%s')
# all_cdr3 = np.loadtxt(ALL_CDR3, fmt='%s')
# 
# #matrix = {'A': pd.DataFrame(index=TRAs), 'B': pd.DataFrame(index=TRBs)}
# matrix = pd.DataFrame(index=all_cdr3, columns=lst)
# 
# F1 = f'tmp_{CHAIN}{CHUNK_I}' # randomly generated filename #'tmp_%s1' %chain
# F2 = ALL_CDR3 # filename of all cdr3 sequences #'tmp_%s2' %chain
# #
# for cdr3_seq in chunk:
#     np.savetxt(F1, np.array([cdr3_seq]), fmt='%s')
#     #
#     cmd = [seq2score, '-blf', BLF, '-blqij', QIJ, '-pa', F1, F2]
#     sim_df = run_cmd(cmd)
#     sim_df = sim_df[:-1]
#     #
#     sim_dct = dict(zip(sim_df.seq2, sim_df.similarity))
#     #sim_dct.update({'': 0}) #Similarity score is 0 when CDR3 seq is ''
#     #
#     matrix[cdr3_seq] = matrix.index.map(sim_dct)
# matrix.to_csv(OUTPUT, index=True)
# 
# #####################
# # END OF SCRIPT TWO #
# #####################
# 
# ######################
# #    SCRIPT THREE    #
# ######################
# import numpy as np
# import pandas as pd
# 
# # Combine chunks and chains
# 
# INPUT = snakemake.input.dfs
# OUTPUT = snakemake.output.df
# 
# dfs = list()
# for filename in INPUTS:
#     tmp = pd.read_csv(filename, index_col=1)
#     dfs.append(tmp, axis=1)
#     
# df = pd.concat(dfs)
# df.loc['missing'] = 0 #Similarity score is 0 when CDR3 seq is ''
# df['missing'] = 0
# 
# assert sum(df.index.duplicated()) == 0
# assert sum(df.index.duplicated()) == 0
# 
# df.to_csv(OUTPUT, index=True)
# 