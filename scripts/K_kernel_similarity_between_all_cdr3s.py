#!/usr/bin/env python

# import numpy as np
# import pandas as pd
# import subprocess
# import os
# 
# def cdr3_lst_converter(x):
#     #define format of datetime
#     return x.replace("[","").replace("]","").replace("'","").split(" ")
# 
# converters={'cdr3_lst_TRA': cdr3_lst_converter, 'cdr3_lst_TRB': cdr3_lst_converter}
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
# IN_FILE = snakemake.params.df #'notebooks/all_contig_annotations.1os2os.csv'
# #'experiments/exp13/run1/tcr/cellranger_all/cellranger_all/outs/multi/vdj_t/all_contig_annotations.csv' #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_TCR/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv' #
# 
# OUTPUTS = snakemake.output.all_cdr3
# 
# # Load data
# df = pd.read_csv(IN_FILE, usecols=['chain','cdr3']) #pd.read_csv(IN_FILE, usecols=['cdr3_lst_TRA', 'cdr3_lst_TRB'], converters=converters, low_memory=False) #
# 
# # OBS! The CDR3 sequences CANNOT contain characters other than the 20 AA letters!
# TRAs = np.unique([seq.strip().strip('*') for seq in df[df.chain == 'TRA'].cdr3.dropna() if seq.strip()!=''])
# TRBs = np.unique([seq.strip().strip('*') for seq in df[df.chain == 'TRB'].cdr3.dropna() if seq.strip()!=''])
# 
# for filename in OUTPUTS:
#     np.savetxt(filename, TRAs, fmt='%s')
# #np.savetxt('all_cdr3_a', TRAs, fmt='%s') # permanently save, tmp_A2
# #np.savetxt('all_cdr3_b', TRBs, fmt='%s') # permanently save
# 
# #####################
# # END OF SCRIPT ONE #
# #####################
#         
# # Write chunck files for parallel processing
# for chain, lst in [('A', TRAs), ('B', TRBs)]:
#     # Get chunks
#     n = 50
#     for i in range(0, len(lst), n):
#         np.savetxt(f'{chain}_chunk{i}', lst[i:i + n], fmt='%s') # permanently save
# 
# 

######################
#     SCRIPT TWO     #
######################

# Run in parallel: Chains and chunks

import numpy as np
import pandas as pd
import subprocess
import os

def run_cmd(cmd, input_string=''):
    """Run the cmd with input_string as stdin and return output."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                             stderr=subprocess.PIPE, universal_newlines=True, close_fds=True)
    result = pd.read_csv(p.stdout, sep=" ", names=['seq1', 'seq2', 'similarity'], usecols=[1,2,3], comment='#')
    p.wait()
    return result

# Input
seq2score = snakemake.input.s2s #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel" #
BLF = snakemake.input.blf #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50" #
QIJ = snakemake.input.qij #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij" #

CHAIN = snakemake.params.chain
CHUNK_I = snakemake.params.chunk_i
CHUNK = snakemake.input.chu
ALL_CDR3 = snakemake.input.cdr

OUTPUT = snakemake.output.chunk

DIR = os.path.dirname(os.path.abspath(CHUNK))

chunk = np.loadtxt(CHUNK, dtype='str')
all_cdr3 = np.loadtxt(ALL_CDR3, dtype='str')

#matrix = {'A': pd.DataFrame(index=TRAs), 'B': pd.DataFrame(index=TRBs)}
matrix = pd.DataFrame(index=all_cdr3, columns=chunk)

F1 = os.path.join(DIR, f'{CHUNK_I}.tmp') # randomly generated filename #'tmp_%s1' %chain
F2 = ALL_CDR3 # filename of all cdr3 sequences #'tmp_%s2' %chain
#
for cdr3_seq in chunk:
    np.savetxt(F1, np.array([cdr3_seq]), fmt='%s')
    #
    cmd = [seq2score, '-blf', BLF, '-blqij', QIJ, '-pa', F1, F2]
    sim_df = run_cmd(cmd)
    sim_df = sim_df[:-1]
    #
    sim_dct = dict(zip(sim_df.seq2, sim_df.similarity))
    #sim_dct.update({'': 0}) #Similarity score is 0 when CDR3 seq is ''
    #
    matrix[cdr3_seq] = matrix.index.map(sim_dct)
matrix.T.to_csv(OUTPUT, index=False, header=False)

#####################
# END OF SCRIPT TWO #
#####################

######################
#    SCRIPT THREE    #
######################
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
#     dfs.append(tmp)
#     
# df = pd.concat(dfs, axis=1)
# df.loc['missing'] = 0 #Similarity score is 0 when CDR3 seq is ''
# df['missing'] = 0
# 
# assert sum(df.index.duplicated()) == 0
# assert sum(df.index.duplicated()) == 0
# 
# df.to_csv(OUTPUT, index=True)
