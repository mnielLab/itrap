#!/usr/bin/env python

import numpy as np
import pandas as pd
import subprocess
import os

def cdr3_lst_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace("'","").split(" ")

converters={'cdr3_lst_TRA': cdr3_lst_converter, 'cdr3_lst_TRB': cdr3_lst_converter}

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
#TMP = snakemake.input.tmp # Hack

IN_FILE = snakemake.params.df #'notebooks/all_contig_annotations.1os2os.csv'
#'experiments/exp13/run1/tcr/cellranger_all/cellranger_all/outs/multi/vdj_t/all_contig_annotations.csv' #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_TCR/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv' #

OUT_FILE_A = snakemake.output.A #'cdr3_a.csv' #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_CAT_IONTORRENT_KMA_AKB/similarity_assessment/cdr3/kernel_similarity_scores/cdr3_a.csv' #
OUT_FILE_B = snakemake.output.B #'cdr3_b.csv' #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_CAT_IONTORRENT_KMA_AKB/similarity_assessment/cdr3/kernel_similarity_scores/cdr3_b.csv' #

# Load data
df = pd.read_csv(IN_FILE, usecols=['chain','cdr3']) #pd.read_csv(IN_FILE, usecols=['cdr3_lst_TRA', 'cdr3_lst_TRB'], converters=converters, low_memory=False) #

# OBS! The CDR3 sequences CANNOT contain characters other than the 20 AA letters!
TRAs = np.unique([seq.strip().strip('*') for seq in df[df.chain == 'TRA'].cdr3.dropna() if seq.strip()!=''])
TRBs = np.unique([seq.strip().strip('*') for seq in df[df.chain == 'TRB'].cdr3.dropna() if seq.strip()!=''])

matrix = {'A': pd.DataFrame(index=TRAs), 'B': pd.DataFrame(index=TRBs)}

np.savetxt('tmp_A2', TRAs, fmt='%s')
np.savetxt('tmp_B2', TRBs, fmt='%s')

for chain, lst in [('A', TRAs), ('B', TRBs)]:
    F1 = 'tmp_%s1' %chain
    F2 = 'tmp_%s2' %chain
    #
    for cdr3_seq in lst:
        np.savetxt(F1, np.array([cdr3_seq]), fmt='%s')
        #
        cmd = [seq2score, '-blf', BLF, '-blqij', QIJ, '-pa', F1, F2]
        sim_df = run_cmd(cmd)
        sim_df = sim_df[:-1]
        #
        sim_dct = dict(zip(sim_df.seq2, sim_df.similarity))
        #sim_dct.update({'': 0}) #Similarity score is 0 when CDR3 seq is ''
        #
        matrix[chain][cdr3_seq] = matrix[chain].index.map(sim_dct)
    matrix[chain].loc['missing'] = 0 #Similarity score is 0 when CDR3 seq is ''
    matrix[chain]['missing'] = 0

assert sum(matrix['A'].index.duplicated()) == 0
assert sum(matrix['B'].index.duplicated()) == 0
assert sum(matrix['A'].columns.duplicated()) == 0
assert sum(matrix['B'].columns.duplicated()) == 0

matrix['A'].to_csv(OUT_FILE_A, index=True)
matrix['B'].to_csv(OUT_FILE_B, index=True)

os.remove("tmp_A1")
os.remove("tmp_A2")
os.remove("tmp_B1")
os.remove("tmp_B2")
