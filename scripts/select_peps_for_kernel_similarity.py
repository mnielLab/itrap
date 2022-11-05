#!/usr/bin/env python

import numpy as np
import pandas as pd
import subprocess

from D_plot_specificity_matrix_utils import epitope_sorter_index

# ARGS
seq2score="/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel"
BLF="/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50"
QIJ="/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij"

#F1 = "/home/tuba/herpov/tcr-pmhc-sc-project/notebooks/lstS.txt"
#F2 = "/home/tuba/herpov/tcr-pmhc-sc-project/notebooks/lstM.txt"
#
#lstS = pd.DataFrame(['CAVRSAYSGAG'])
#lstM = pd.DataFrame(['CAARLIQGAQKLVF', 'CAGPSYNTDKLIF', 'CAMPNSGGYQKVTF', 'CAMNRDDKIIF', 'CAVRSAYSGAGSYQLTF'])

EXP = "exp3"
PLATFORM = "IONTORRENT"
MAPPING = 'KMA' # BLAST
BARCODE_SYSTEM = 'AKB'

SRT_FIL = ("/home/tuba/herpov/tcr-pmhc-sc-project/data/" +
           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
           "/specificity_matrix/peptide_per_clonotype_by_gem_size/unique_tcrs/v1.b1.t1.ecs_False.ess_False.lst")

GEM_FIL = ("/home/tuba/herpov/tcr-pmhc-sc-project/data/" +
           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
           "/specificity_matrix/peptide_per_clonotype_by_gem_size/unique_gems/v1.b1.t1.ecs_False.ess_False.lst")

IN_FILE = ("/home/tuba/herpov/tcr-pmhc-sc-project/data/" +
           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
           "/tables/tcr_barcode.augmented.csv")

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
df = pd.read_csv(IN_FILE, usecols=['gem', 'clonotype', 'num_clonotype', 'cdr3_TRA', 'cdr3_TRB', 'epitope', 'peptide', 'read_counts_mhc', 'binding_concordance'])
sorted_tcrs = np.loadtxt(SRT_FIL, dtype='int8')
unique_gems = np.loadtxt(GEM_FIL, dtype='U20')

# MAIN
df = df[df.gem.isin(unique_gems)]
df['rank'] = sorter_index(sorted_tcrs, df.num_clonotype)
df.sort_values(by=['rank'], inplace=True)
df.replace({'cdr3_TRA': "", 'cdr3_TRB': ""}, np.nan, inplace=True)
df.dropna(subset=(['cdr3_TRA', 'cdr3_TRB']), inplace=True)
df['peptide'] = df.peptide.str.split("_", expand=True)[0].str.strip()
#df['epitope_rank'] = epitope_sorter_index(df)
#df.sort_values(by=['epitope_rank', 'clonotype'], inplace=True)
df.reset_index(drop=True, inplace=True)

matrix = df[['gem', 'peptide']].copy()

np.savetxt('tmp_P2', df.peptide.values, fmt='%s')

for index, row in df.iterrows():
    pep = row.peptide
    
    np.savetxt('tmp_P1', np.array([row.peptide]), fmt='%s')
    
    cmd = [seq2score, '-blf', BLF, '-blqij', QIJ, '-pa', 'tmp_P1', 'tmp_P2']
    sim_df = run_cmd(cmd)
    sim_df = sim_df[:-1]

    assert sum(sim_df.seq2.values == df.peptide.values) == df.shape[0], "The sequences are not listed similarly"
    #assert sim_df.shape[0] == matrix[chain].shape[0], "The sizes do not match"
    assert sim_df.seq1.unique()[0] == row.peptide, "The tested sequence does not match the CDR3 sequence"

    matrix = pd.concat([matrix, sim_df[['similarity']]], axis=1, sort=False) #

    assert df.shape[0] == matrix.shape[0], "The sizes do not match"       
    matrix.rename(columns={'similarity': row.gem}, inplace=True)

matrix.set_index(['gem'], inplace=True)

df.set_index('gem', inplace=True)
assert sum(matrix.index == df.index) == df.shape[0], "The indexes are not consistent"
matrix = pd.merge(df[['clonotype', 'epitope', 'peptide', 'read_counts_mhc', 'binding_concordance']], matrix, on='peptide', left_index=True, right_index=True)
assert df.shape[0] == matrix.shape[0], "The sizes do not match"

matrix.to_csv(OUT_DIR + 'peptide.csv', index=True) # GEMs are indexes, thus DO write index to file


