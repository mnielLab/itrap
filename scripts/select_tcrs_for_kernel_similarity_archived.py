#!/usr/bin/env python

import numpy as np
import pandas as pd
import subprocess

#from D_plot_specificity_matrix_utils import epitope_sorter_index

# ARGS
seq2score="/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel"
BLF="/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50"
QIJ="/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij"

#F1 = "/home/tuba/herpov/tcr-pmhc-sc-project/notebooks/lstS.txt"
#F2 = "/home/tuba/herpov/tcr-pmhc-sc-project/notebooks/lstM.txt"
#
#lstS = pd.DataFrame(['CAVRSAYSGAG'])
#lstM = pd.DataFrame(['CAARLIQGAQKLVF', 'CAGPSYNTDKLIF', 'CAMPNSGGYQKVTF', 'CAMNRDDKIIF', 'CAVRSAYSGAGSYQLTF'])

EXP = "exp6"
PLATFORM = "IONTORRENT"
MAPPING = 'KMA' # BLAST
BARCODE_SYSTEM = 'AKB'

SRT_FIL = ("/home/tuba/herpov/tcr-pmhc-sc-project/data/" +
           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
           "/specificity_matrix/peptide_per_clonotype_by_gem_size/unique_tcrs/b1.t1.ecs_False.ess_False.lst") #v1.

GEM_FIL = ("/home/tuba/herpov/tcr-pmhc-sc-project/data/" +
           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
           "/specificity_matrix/peptide_per_clonotype_by_gem_size/unique_gems/b1.t1.ecs_False.ess_False.lst") #v1.

IN_FILE = ("/home/tuba/herpov/tcr-pmhc-sc-project/data/" +
           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
           "/tables/tcr_barcode.cleaned.csv")

OUT_DIR = ("/home/tuba/herpov/tcr-pmhc-sc-project/data/" +
           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM + "/similarity_assessment/")


def run_cmd(cmd, input_string=''):
    """Run the cmd with input_string as stdin and return output."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                             stderr=subprocess.PIPE, universal_newlines=True, close_fds=True)
    result = pd.read_csv(p.stdout, sep=" ", names=['seq1', 'seq2', 'similarity'], usecols=[1,2,3], comment='#')
    p.wait()
    return result

def sorter_index(sorter, var):
    sorterIndex = dict(zip(sorter,range(len(sorter))))
    return var.map(sorterIndex)

# Load data
df = pd.read_csv(IN_FILE, usecols=['gem', 'ct', 'clonotype', 'num_clonotype', 'cdr3_TRA', 'cdr3_TRB', 'epitope', 'peptide', 'umi_count_mhc', 'umi_count_tcr', 'binding_concordance'])
sorted_tcrs = np.loadtxt(SRT_FIL, dtype='int8')
unique_gems = np.loadtxt(GEM_FIL, dtype='U20')

# MAIN
df = df[df.gem.isin(unique_gems)]
df['rank'] = sorter_index(sorted_tcrs, df.num_clonotype)
df.sort_values(by=['rank'], inplace=True)
df.replace({'cdr3_TRA': "", 'cdr3_TRB': ""}, np.nan, inplace=True)
df.dropna(subset=(['cdr3_TRA', 'cdr3_TRB']), inplace=True)
#df['epitope_rank'] = epitope_sorter_index(df)
#df.sort_values(by=['epitope_rank', 'clonotype'], inplace=True)
df.reset_index(drop=True, inplace=True)

matrix = {'A': df[['gem', 'cdr3_TRA']].copy(), 'B': df[['gem', 'cdr3_TRB']].copy()}

np.savetxt('tmp_A2', df.cdr3_TRA.values, fmt='%s')
np.savetxt('tmp_B2', df.cdr3_TRB.values, fmt='%s')

for index, row in df.iterrows():
    tcr = row.clonotype
    
    # assert that each clonotype only contains a single combination of TRA and TRB
    # This is false, but it ought to be that way?
    #assert df[df.clonotype==tcr].groupby(['cdr3_TRA', 'cdr3_TRB']).size().values.shape == (1,), "%s was annotated with multiple cdr3 chains" %tcr
    if not df[df.clonotype==tcr].groupby(['cdr3_TRA', 'cdr3_TRB']).size().values.shape == (1,): print(tcr)

    for chain in ['A', 'B']:
        F1 = 'tmp_%s1' %chain
        F2 = 'tmp_%s2' %chain
        cdr3_name = 'cdr3_TR%s' %chain
        
        np.savetxt(F1, np.array([row[cdr3_name]]), fmt='%s') #np.savetxt('tmp_A1', np.array([row.cdr3_TRA]), fmt='%s')
        
        cmd = [seq2score, '-blf', BLF, '-blqij', QIJ, '-pa', F1, F2]
        sim_df = run_cmd(cmd)
        sim_df = sim_df[:-1]

        assert sum(sim_df.seq2.values == df[cdr3_name].values) == df.shape[0], "The sequences are not listed similarly"
        #assert sim_df.shape[0] == matrix[chain].shape[0], "The sizes do not match"
        assert sim_df.seq1.unique()[0] == row[cdr3_name], "The tested sequence does not match the CDR3 sequence"

        matrix[chain] = pd.concat([matrix[chain], sim_df[['similarity']]], axis=1, sort=False) #

        assert df.shape[0] == matrix[chain].shape[0], "The sizes do not match"       
        matrix[chain].rename(columns={'similarity': row.gem}, inplace=True)

matrix['A'].set_index(['gem'], inplace=True)
matrix['B'].set_index(['gem'], inplace=True)

matrix['A'].rename(columns={'cdr3_TRA': 'cdr3'}, inplace=True)
matrix['B'].rename(columns={'cdr3_TRB': 'cdr3'}, inplace=True)

matrix_C = matrix['A'].add(matrix['B'], fill_value=0)

df.set_index('gem', inplace=True)
assert sum(matrix_C.index == df.index) == df.shape[0], "The indexes are not consistent"
matrix_C = pd.merge(df[['ct', 'clonotype', 'num_clonotype', 'epitope', 'peptide', 'umi_count_mhc', 'umi_count_tcr', 'binding_concordance']], matrix_C, left_index=True, right_index=True)
assert df.shape[0] == matrix_C.shape[0], "The sizes do not match"

matrix_C.to_csv(OUT_DIR + 'cdr3.csv', index=True) # GEMs are indexes, thus DO write index to file
matrix['A'].to_csv(OUT_DIR + 'cdr3_a.csv', index=True)
matrix['B'].to_csv(OUT_DIR + 'cdr3_b.csv', index=True)


