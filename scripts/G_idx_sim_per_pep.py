#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from ast import literal_eval
import re
from scipy import stats
from random import sample
#import scipy.stats as stats
import random
import re
import statistics
from scipy.stats import wilcoxon
from sklearn.metrics import roc_curve, auc
import os

#plt.style.use('ggplot')

def HLA_cd8_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace(",", "").replace("'","").split(" ")

def cdr3_lst_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace("'","").split(" ")

def epitope_converter(x):
    #define format of datetime
    return [y for y in x.replace("[","").replace("]","").replace("\n","").split("'") if (y != '') & (y != ' ')]

def peptide_hla_converter(x):
    return re.findall("\w+\s{1}\w{1}\d+", x.replace("[","").replace("]","").replace("\n","").replace("'",""))

def literal_converter(val):
    # replace NaN with '' and perform literal eval on the rest
    return [] if val == '' else literal_eval(val)

converters = {'peptide_HLA_lst': peptide_hla_converter,
              'umi_count_lst_mhc': literal_eval,
              'umi_count_lst_TRA': literal_converter,'umi_count_lst_TRB': literal_converter,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'cdr3_lst_TRB': cdr3_lst_converter,
              'HLA_lst_mhc': cdr3_lst_converter,'HLA_cd8': HLA_cd8_converter} #


import sys  
sys.path.insert(0, '../scripts')

from D_plot_specificity_matrix_utils import calc_binding_concordance


###########################################################
#                        Functions                        #
###########################################################
def keep_max_conc(df):
    df = calc_binding_concordance(df.copy(), 'ct')
    assert df.binding_concordance.isna().any() == False
    
    dct = df.groupby(['ct']).binding_concordance.max()
    df['max_conc'] = df.ct.map(df.groupby(['ct']).binding_concordance.max())
    return df.binding_concordance == df.max_conc

def get_filtered_data(df, opt_thr):
    # Convert series of opt_thr to an exacutable and execute
    selection = eval(' & '.join(['(df.%s >= %d)'%(k,v) for k,v in opt_thr.iterrows()]))
    return selection

def get_multiplets(df):
    #tmp = df[idx1 & idx2]
    dct = df.groupby(['ct','peptide_HLA']).gem.count() > 1
    idx = df.set_index(['ct','peptide_HLA']).index.map(dct)
    return idx.fillna(False)

def subset_similarity_scores(i, tmp_filename): # append intra or inter to filename
    #if len(i) > 1000:
    #    i = np.sort(np.random.choice(i, 1000, replace=False))
    np.savetxt(tmp_filename, i+1, fmt='%d')
    #cmd = f"sed -n '{'p;'.join(map(str,i+1))}p' {filename}" #sed -n '1p;3p' # OBS! sed is not zero-indexed!
    #
    #f = open(tmp_filename, "w")
    #f.write(cmd)
    #f.close()
    
def get_index(ci,cj,chain_idx):
    if type(ci) is str:
        i = np.arange(len(chain_idx))[chain_idx == ci]
        assert len(i) == 1, i
    else:
        i = np.nonzero(ci[:,None] == chain_idx)[1]
        assert len(ci) == len(i), ci
        
    if type(cj) is str:
        j = np.arange(len(chain_idx))[chain_idx == cj]
        assert len(j) == 1, j
    else:
        j = np.nonzero(cj[:,None] == chain_idx)[1]
        assert len(cj) == len(j), cj
        
    return i,j

def split_similarity_scores(ai, aj, bi, bj): # input intra or inter
    """
    plateau : intra
    ai and bi are lists of strings
    aj and bj are lists of strings
    
    plateau : inter
    ai and bi are lists of strings
    aj and bj are strings
    
    indexes: find the indexes anything by the elements of ai and bi
    columns: find the index of the single column containing aj/bj
    """
    assert len(ai) == len(bi)
    
    i_tra, j_tra = get_index(ai, aj, idx_tra)   
    i_trb, j_trb = get_index(bi, bj, idx_trb)
    
    p = '_'.join(peptide.split())
    
    tmp_a = os.path.join(SIM_OUTPUT, f'a.{p}.txt')
    tmp_b = os.path.join(SIM_OUTPUT, f'b.{p}.txt')
    
    subset_similarity_scores(np.unique(i_tra), tmp_a)
    subset_similarity_scores(np.unique(i_trb), tmp_b)

def notnan(x):
    return x == x

##########################################################
#                         Inputs                         #
##########################################################
#CAT_DIR = '/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run2/cat/'
#CAT_DIR = '../experiments/exp13/run1_archive/cat/'

VALID = snakemake.input.valid 
#SIM_TRA = snakemake.input.sim_tra #CAT_DIR + "/similarity_assessment/cdr3_a.csv"
IDX_TRA = snakemake.input.idx_tra
#SIM_TRB = snakemake.input.sim_trb #CAT_DIR + "/similarity_assessment/cdr3_b.csv"
IDX_TRB = snakemake.input.idx_trb

##########################################################
#                         Output                         #
##########################################################
SIM_OUTPUT = snakemake.output[0]

##########################################################
#                          Load                          #
##########################################################
print('loading data')
df = pd.read_csv(VALID, converters=converters, low_memory=False)

idx_tra = np.loadtxt(IDX_TRA, dtype='str')
idx_trb = np.loadtxt(IDX_TRB, dtype='str')

##########################################################
#                          Prep                          #
##########################################################
# For the similarity computation we don't want to include clonotypes that are missing a chain
df.dropna(subset=['cdr3_TRA','cdr3_TRB'], inplace=True)
# Remove GEMs with weird CDR3 annotations
expected_tra = df.cdr3_TRA.str.contains('^[GPAVLIMCFYWHKRQNEDST]+$', regex=True, na=False)
expected_trb = df.cdr3_TRB.str.contains('^[GPAVLIMCFYWHKRQNEDST]+$', regex=True, na=False)
df = df[expected_tra & expected_trb] # This step should be done initially when cleaning the TCR sequences.

df.replace('0','', inplace=True) # why?
df.sort_values(by='epitope_rank', inplace=True)

##########################################################
#                   Index per peptide                  #
##########################################################
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
   
sdf = df[keep_max_conc(df)]
grp = sdf.groupby('peptide_HLA', sort=False)
for i, (peptide, group) in enumerate(grp):
    group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB'], inplace=True)
    if len(group) == 1:
        continue
    if len(group) > 1000:
        group = group.sample(n=1000, replace=False, random_state=1)
        

    inter_chains = sdf.loc[sdf.peptide_HLA != peptide, ['cdr3_TRA', 'cdr3_TRB']] 

    # OBS! Maybe get the pairs directly instead of zipping them?
    # Make sure you get true pairs and not just random pairs?!
    #cdr3_TRAs = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRA.values
    #cdr3_TRBs = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRB.values
    cdr3_TRAs = group.cdr3_TRA.values
    cdr3_TRBs = group.cdr3_TRB.values

    assert len(cdr3_TRAs) == len(cdr3_TRBs) == len(group.loc[:,['cdr3_TRA','cdr3_TRB']]) #.drop_duplicates()
    
    print(i, peptide, len(cdr3_TRAs))

    split_similarity_scores(cdr3_TRAs, cdr3_TRAs, cdr3_TRBs, cdr3_TRBs)