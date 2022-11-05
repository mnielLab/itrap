#!/usr/bin/env python

import numpy as np
import pandas as pd
import subprocess

#from D_plot_specificity_matrix_utils import epitope_sorter_index

# ARGS
seq2score=snakemake.input.s2s #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel"
BLF=snakemake.input.blf #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50"
QIJ=snakemake.input.qij #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij"

IN_FILE = snakemake.input.df

OUT_FILE_C = snakemake.output.C
OUT_FILE_A = snakemake.output.A
OUT_FILE_B = snakemake.output.B


def run_cmd(cmd, input_string=''):
    """Run the cmd with input_string as stdin and return output."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                             stderr=subprocess.PIPE, universal_newlines=True, close_fds=True)
    result = pd.read_csv(p.stdout, sep=" ", names=['seq1', 'seq2', 'similarity'], usecols=[1,2,3], comment='#')
    p.wait()
    return result

# Load data
df = pd.read_csv(IN_FILE, usecols=['gem', 'ct', 'clonotype', 'num_clonotype', 'cdr3_TRA', 'cdr3_TRB', 'epitope', 'peptide', 'umi_count_mhc', 'umi_count_tcr', 'binding_concordance'], low_memory=False)

## MAIN
#df = df[df.gem.isin(unique_gems)]
df.fillna({'cdr3_TRA':'', 'cdr3_TRB':''}, inplace=True)
#df['rank'] = sorter_index(sorted_tcrs, df.num_clonotype)
#df.sort_values(by=['rank'], inplace=True)
#df.replace({'cdr3_TRA': "", 'cdr3_TRB': ""}, np.nan, inplace=True)
#df.dropna(subset=(['cdr3_TRA', 'cdr3_TRB']), inplace=True)
##df['epitope_rank'] = epitope_sorter_index(df)
##df.sort_values(by=['epitope_rank', 'clonotype'], inplace=True)
df.reset_index(drop=True, inplace=True)

matrix = {'A': df[['gem', 'cdr3_TRA']].copy(), 'B': df[['gem', 'cdr3_TRB']].copy()}

np.savetxt('tmp_A2', np.array(list(filter(None, df.cdr3_TRA.unique()))), fmt='%s') #values
np.savetxt('tmp_B2', np.array(list(filter(None, df.cdr3_TRB.unique()))), fmt='%s') #values

for index, row in df.iterrows():
    tcr = row.clonotype
    
    # assert that each clonotype only contains a single combination of TRA and TRB
    # This is false, but it ought to be that way?
    #assert df[df.clonotype==tcr].groupby(['cdr3_TRA', 'cdr3_TRB']).size().values.shape == (1,), "%s was annotated with multiple cdr3 chains" %tcr
    if not df[df.clonotype==tcr].groupby(['cdr3_TRA', 'cdr3_TRB']).size().values.shape[0] == 1: print(tcr)

    for chain in ['A', 'B']:
        F1 = 'tmp_%s1' %chain
        F2 = 'tmp_%s2' %chain
        cdr3_name = 'cdr3_TR%s' %chain

        if row[cdr3_name] == '':
            matrix[chain][row.gem] = 0 #Similarity score is 0 when CDR3 seq is ''
            break
        
        np.savetxt(F1, np.array([row[cdr3_name]]), fmt='%s') #np.savetxt('tmp_A1', np.array([row.cdr3_TRA]), fmt='%s')
        
        cmd = [seq2score, '-blf', BLF, '-blqij', QIJ, '-pa', F1, F2]
        sim_df = run_cmd(cmd)
        sim_df = sim_df[:-1]
        
        sim_dct = dict(zip(sim_df.seq2, sim_df.similarity))
        sim_dct.update({'': 0}) #Similarity score is 0 when CDR3 seq is ''
        #print(sim_dct)

        matrix[chain][row.gem] = matrix[chain][cdr3_name].map(sim_dct)

        #assert sum(sim_df.seq2.values == df[cdr3_name].values) == df.shape[0], "The sequences are not listed similarly" #
        ##assert sim_df.shape[0] == matrix[chain].shape[0], "The sizes do not match"
        #assert sim_df.seq1.unique()[0] == row[cdr3_name], "The tested sequence does not match the CDR3 sequence"

        #matrix[chain] = pd.concat([matrix[chain], sim_df[['similarity']]], axis=1, sort=False) #

        #assert df.shape[0] == matrix[chain].shape[0], "The sizes do not match"       
        #matrix[chain].rename(columns={'similarity': row.gem}, inplace=True)

matrix['A'].set_index(['gem'], inplace=True)
matrix['B'].set_index(['gem'], inplace=True)

matrix['A'].rename(columns={'cdr3_TRA': 'cdr3'}, inplace=True)
matrix['B'].rename(columns={'cdr3_TRB': 'cdr3'}, inplace=True)

matrix_C = matrix['A'].add(matrix['B'], fill_value=0)

cols = matrix_C.columns.tolist() #matrix_C.columns.tolist().remove('cdr3')
cols.remove('cdr3')
df.set_index('gem', inplace=True)
assert sum(matrix_C.index == df.index) == df.shape[0], "The indexes are not consistent"
matrix_C = pd.merge(df[['ct', 'clonotype', 'num_clonotype', 'epitope', 'peptide', 'umi_count_mhc', 'umi_count_tcr', 'binding_concordance', 'cdr3_TRA', 'cdr3_TRB']],
                    matrix_C[['cdr3'] + cols],
                    left_index=True, right_index=True)
assert df.shape[0] == matrix_C.shape[0], "The sizes do not match"

matrix_C.to_csv(OUT_FILE_C, index=True) # GEMs are indexes, thus DO write index to file
matrix['A'].to_csv(OUT_FILE_A, index=True)
matrix['B'].to_csv(OUT_FILE_B, index=True)


