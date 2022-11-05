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

def get_intra_chains(grp, pep):
    # OBS! Maybe get the pairs directly instead of zipping them?
    # Make sure you get true pairs and not just random pairs?!
    cdr3_TRAs = grp.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRA.values
    cdr3_TRBs = grp.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRB.values
    
    return cdr3_TRAs, cdr3_TRBs

def get_unique_entries(df):
    return np.where((df.cdr3_TRA==a) & (df.cdr3_TRB==b), False, True)

def get_sample_size():
    inter_entries = get_unique_entries(inter_chains)
    inter_indexes = inter_chains[inter_entries].index.to_list()
    #print(sum(inter_entries), len(inter_indexes), sum(get_unique_entries(group)))
    return min(sum(get_unique_entries(group)), len(inter_indexes))
    
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

def add_similarity_scores(ai,a, bi,b, plateau):
    if plateau == 'inter':
        #mat_a, mat_b = load_similarity_scores(ai,a, bi,b, plateau)
        
        ai, a = get_index(ai, a, idx_tra)
        bi, b = get_index(bi, b, idx_trb)
        
        #mat_a = mat_a.loc[ai, a].reset_index(drop=True).T.reset_index(drop=True).T
        #mat_b = mat_b.loc[bi, b].reset_index(drop=True).T.reset_index(drop=True).T
        
        mat_a = sim_tra.loc[ai, a].reset_index(drop=True).T.reset_index(drop=True).T
        mat_b = sim_trb.loc[bi, b].reset_index(drop=True).T.reset_index(drop=True).T
        
        #sim_tra, sim_trb = load_similarity_scores(ai,a, bi,b, plateau)
        #mat_a = sim_tra.reset_index(drop=True).T.reset_index(drop=True).T
        #mat_b = sim_trb.reset_index(drop=True).T.reset_index(drop=True).T
        
        assert len(mat_a) == len(mat_b) == len(ai) == len(bi)
        assert len(a) == len(b) == mat_a.shape[1] == mat_b.shape[1]
        
    elif plateau == 'intra':
        ai, a = get_index(ai, a, idx_tra)
        bi, b = get_index(bi, b, idx_trb)
        
        mat_a = sim_tra.loc[ai, a].reset_index(drop=True).T.reset_index(drop=True).T
        mat_b = sim_trb.loc[bi, b].reset_index(drop=True).T.reset_index(drop=True).T
        
        assert len(mat_a) == len(mat_b)
        assert len(ai) == len(bi)
        
        
    return mat_a.add(mat_b)

def avg_similarity_scores(ai,bi):
    mat_a = sim_tra.loc[ai, a].reset_index(drop=True).T.reset_index(drop=True).T.to_frame().T
    mat_b = sim_trb.loc[bi, b].reset_index(drop=True).T.reset_index(drop=True).T.to_frame().T
    return pd.concat([mat_a, mat_b], ignore_index=True).mean() #.to_frame() is important?

def get_intra_similarity(cdr3_TRAs, cdr3_TRBs, method='sum'):
    unique_entries = get_unique_entries(group)
    unique_entry_indexes = group[unique_entries].index.to_list()

    sample_size = get_sample_size() 
    # Am I cheating myself by having the chance of not sampling the most similar chains?
    # There can be quite substantial size differences: 100x
    
    sampled_idxs = random.sample(unique_entry_indexes, sample_size)
    
    intra_a = group.loc[sampled_idxs, 'cdr3_TRA'].values
    intra_b = group.loc[sampled_idxs, 'cdr3_TRB'].values
    
    if method == 'sum':
        #combined_similarity = add_similarity_scores(intra_a,a, intra_b,b, 'intra')
        combined_similarity = add_similarity_scores(a,intra_a, b,intra_b, 'intra') #, reversed index and columns!
    else:
        combined_similarity = avg_similarity_scores(intra_a, intra_b)
    return {'score': combined_similarity.max(axis=1), #added axis=1
            'fraction': sum(combined_similarity > 1.8)/len(combined_similarity)}

def get_inter_similarity(cdr3_TRAs, cdr3_TRBs, method='sum'):
    # OBS! make sure the size to sample from matches the number og unique entries intra_similarity! 
    unique_entries = get_unique_entries(inter_chains)
    unique_entry_indexes = inter_chains[unique_entries].index.to_list()

    sample_size = get_sample_size() # Am I cheating myself by having the chance of not sampling the most similar chains?

    sampled_idxs = random.sample(unique_entry_indexes, sample_size)
    inter_a = inter_chains.loc[sampled_idxs, 'cdr3_TRA'].values
    inter_b = inter_chains.loc[sampled_idxs, 'cdr3_TRB'].values
    
    if method == 'sum':
        #combined_similarity = add_similarity_scores(inter_a,a, inter_b,b, 'inter')
        combined_similarity = add_similarity_scores(a,inter_a, b,inter_b, 'inter')
    else:
        combined_similarity = avg_similarity_scores(inter_a, inter_b)
    return {'score': combined_similarity.max(axis=1), # added axis=1
            'fraction': sum(combined_similarity > 1.8)/len(combined_similarity)}

def add_data(scores, peptide, hue, filtering, random_sample):
    n = len(scores)
    pep_lst = [peptide] * n
    hue_lst = [hue] * n
    flt_lst = [filtering] * n
    rnd_lst = [random_sample] * n
    
    #g = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB'])
    
    #cdr3s = list(zip(g.cdr3_TRA.values, g.cdr3_TRB.values))
    cdr3s = list(zip(cdr3_TRAs, cdr3_TRBs))
    cts = g.ct.values
    
    assert len(cdr3s) == n
    assert len(cts) == n
    
    tmp = pd.DataFrame(np.array([scores, pep_lst, hue_lst, flt_lst, rnd_lst, cdr3s, cts]).T,
                       columns=['score','peptide','plateau','filtering', 'rnd_sample','cdr3s','ct'])
    tmp.score = tmp.score.astype(float)
    tmp.set_index(['plateau','rnd_sample','filtering','peptide','ct'], inplace=True)
    #tmp.index = pd.MultiIndex.from_product([plateaus, random_sampling, filterings, peptides, clonotypes],
    #                                       names=['plateau', 'rnd_sample','filtering','peptide','ct'])
    return tmp

def notnan(x):
    return x == x

##########################################################
#                         Inputs                         #
##########################################################
#CAT_DIR = '/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run2/cat/'
#CAT_DIR = '../experiments/exp13/run1_archive/cat/'

VALID = snakemake.input.valid
THRESHOLD = snakemake.input.opt_thr
GEX = snakemake.input.gex #'../experiments/exp13/run1_archive/tcr/usable_gems.txt'
#TCR = '../experiments/exp13/run1/tcr/cellranger_tot/outs/multi/vdj_t/all_contig_annotations.csv'
#TCR_ARC = '../experiments/exp13/run1_archive/tcr/cellranger_tot/outs/multi/vdj_t/all_contig_annotations.csv'

SIM_TRA = snakemake.input.sim_tra
SIM_TRB = snakemake.input.sim_trb
IDX_TRA = snakemake.input.idx_tra
IDX_TRB = snakemake.input.idx_trb
IDX_TRA_SUB = snakemake.input.idx_tra_sub
IDX_TRB_SUB = snakemake.input.idx_trb_sub



filter_set = snakemake.params.flt
random_sample = int(snakemake.params.rnd)
peptide = ' '.join(snakemake.params.pep.split('_'))

##########################################################
#                         Output                         #
##########################################################
SIM_OUTPUT = snakemake.output.sim
#AUC_OUTPUT = snakemake.output.auc

##########################################################
#                          Load                          #
##########################################################

#hto = pd.read_csv(HTO, skiprows=1, header=None,
#                  names=['gem','seurat','umi_count_hto','feature_rna','count_hto','feature_hto',
#                         'hto_max_id','hto_sec_id','hto_margin','hto_classification','hto_global_class','hash_id'])
#tcr = pd.read_csv(TCR)
#tcr_dct = tcr.groupby('barcode').is_cell.unique().apply(lambda x: x[0])
#tcr_arc = pd.read_csv(TCR_ARC)
#arc_dct = tcr_arc.groupby('barcode').is_cell.unique().apply(lambda x: x[0])
#tcr_cell = pd.merge(arc_dct,tcr_dct, left_index=True,right_index=True, how='outer', suffixes=['_arc','_gex'])

opt_thr = pd.read_csv(THRESHOLD, index_col=0, header=None, names=['thr']).thr.dropna()
gex = pd.read_csv(GEX, header=None, names=['gem'])
print('loading data')
df = pd.read_csv(VALID, converters=converters, low_memory=False)


##########################################################
#                          Prep                          #
##########################################################
# when filtering for optimal threshold it is important to have values in UMI and delta
df.fillna({'umi_count_mhc':0, 'delta_umi_mhc':0, "umi_count_mhc_rel":0,
           'umi_count_TRA':0, 'delta_umi_TRA':0,
           'umi_count_TRB':0, 'delta_umi_TRB':0}, inplace=True)

# Add extra features
df['gex'] = df.gem.isin(gex.gem)
df.single_barcode_mhc = np.where(df.single_barcode_mhc, 'pMHC singlet','pMHC multiplet')
df['clonotype_multiplet'] = df.ct.map(df.groupby('ct').size() > 1)
df['HLA_match_per_gem'] = df.apply(lambda row: row.HLA_mhc in row.HLA_cd8 if row.HLA_cd8 == row.HLA_cd8 else False, axis=1)

print(df.shape)
# For the similarity computation we don't want to include clonotypes that are missing a chain
df.dropna(subset=['cdr3_TRA','cdr3_TRB'], inplace=True)
print(df.shape)
# Remove GEMs with weird CDR3 annotations
expected_tra = df.cdr3_TRA.str.contains('^[GPAVLIMCFYWHKRQNEDST]+$', regex=True, na=False)
expected_trb = df.cdr3_TRB.str.contains('^[GPAVLIMCFYWHKRQNEDST]+$', regex=True, na=False)
df = df[expected_tra & expected_trb] # This step should be done initially when cleaning the TCR sequences.
print(df.shape)
df.replace('0','', inplace=True) # why?
df.sort_values(by='epitope_rank', inplace=True)

##########################################################
#                   Prep similarity load                 #
##########################################################
#cdr3_TRAs, cdr3_TRBs = get_intra_chains(df[keep_max_conc(df)].groupby('peptide_HLA').get_group(peptide), peptide)

idx_tra = np.loadtxt(IDX_TRA, dtype='str')
idx_trb = np.loadtxt(IDX_TRB, dtype='str')

idx_tra_sub = np.loadtxt(IDX_TRA_SUB, dtype='str').astype(int) - 1 # convert from 1-indexed to 0-indexed
idx_trb_sub = np.loadtxt(IDX_TRB_SUB, dtype='str').astype(int) - 1

assert sum(idx_tra_sub < 0) == 0
assert sum(idx_trb_sub < 0) == 0

print(len(idx_tra_sub),len(idx_trb_sub))

sim_tra = pd.read_csv(SIM_TRA, header=None)
sim_trb = pd.read_csv(SIM_TRB, header=None)

print(len(sim_tra), len(sim_trb))

assert len(idx_tra_sub) == len(sim_tra)
assert len(idx_trb_sub) == len(sim_trb)

#print('getting _indexes')
#i_tra, _ = get_index(cdr3_TRAs, cdr3_TRAs, idx_tra)   
#i_trb, _ = get_index(cdr3_TRBs, cdr3_TRBs, idx_trb)
#
## the number of rows should match the number of unique queries
#assert len(sim_tra) == len(np.unique(i_tra)) == len(np.unique(cdr3_TRAs))
#assert len(sim_trb) == len(np.unique(i_trb)) == len(np.unique(cdr3_TRBs))
## The number of queries from a and b should be identical
#assert len(i_tra) == len(i_trb)
#assert idx_tra == np.unique(i_tra)
#assert idx_trb == np.unique(i_trb)

sim_tra.index = idx_tra_sub #np.unique(i_tra)
sim_trb.index = idx_trb_sub #np.unique(i_trb)

print('loaded sim_tra', sim_tra.shape)
print('loaded sim_trb', sim_trb.shape)



##########################################################
#                         Filters                        #
##########################################################
idx0 = ~df.gem.isna() # Total
idx1 = eval(' & '.join([f'(df.{k} >= {abs(v)})' for k,v in opt_thr.items()])) # optimal threshold
idx2 = df.hto_global_class == 'Singlet' # Hashing singlets
idx3 = df.apply(lambda row: row.HLA_mhc in row.HLA_cd8 if (notnan(row.peptide_HLA) & notnan(row.HLA_cd8)) else False, axis=1) # Matching HLA
idx4 = df['exclude_single-chain_TCRs'] # Complete TCRs
idx5 = get_multiplets(df) # Only specificities in multiplets
try:
    idx6 = df.cell_flag # is_cell
except AttributeError:
    idx6 = df.gem.isna() # All false
#idx7 = df.is_cell_gex
idx8 = df.gex # Gene expression data
idxX = ~df.peptide_HLA.isin(['SLEGGGLGY A0101', 'STEGGGLAY A0101', 'ALIAPVHAV A0201',
                             'AYSSAGASI A2402', 'GPAESAAGL B0702', 'AAKGRGAAL B0801']) # 10x negatives
try:
    idx9 = ~df.peptide_HLA_icon.isna() #~df.gem_icon.isna()
    idx10 = df.apply(lambda row: row.HLA_mhc_icon in row.HLA_cd8, axis=1)
except AttributeError:
    idx9 = df.gem.isna()
    idx10 = df.gem.isna()

is_cell_gex = ' (GEX)' if any(idx8) else ''

print('testing filters')
if filter_set == 'indv':
    # Showing individual effects of filtering
    filterings = [idx0,
                  idxX,
              idx1,
              idx3,
              idx2,
              idx4,
              idx5,
              idx6,
              #idx7,
              idx8, idx9, idx10]
    labels = ['total','no negatives','optimal threshold',
          'matching HLA',
          'hashing singlets',
          'complete TCRs',
          'specificity multiplets',
          'is cell%s' %is_cell_gex,
          'is viable cell', 'icon','icon matching HLA']
    palette = ['grey','#ccad60','yellow','#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','black','#d6b4fc','violet']
    
    flt_to_remove = list()
    for i, flt in enumerate(filterings):
        print(labels[i], sum(flt))
        if sum(flt) == 0:
            print('removed:', labels[i])
            flt_to_remove.append(i)
            
    for i in flt_to_remove[::-1]:
        del filterings[i]
        del labels[i]
        del palette[i]

elif filter_set == 'comb':
    # Showing combined effects in the same order
    #filterings = [idx0,
    #              idxX,
    #              (idx1),
    #              (idx1 & idx3),
    #              (idx1 & idx2 & idx3),
    #              (idx1 & idx2 & idx3 & idx4),
    #              (idx1 & idx2 & idx3 & idx4 & idx5),
    #              (idx1 & idx2 & idx3 & idx4 & idx5 & idx6),
    #              #(idx1 & idx2 & idx3 & idx4 & idx5 & idx7),
    #              (idx1 & idx2 & idx3 & idx4 & idx5 & idx8),
    #              idx9, (idx9 & idx10)] # ICON
    labels = ['total','no negatives','optimal threshold',
              'matching HLA',
              'hashing singlets',
              'complete TCRs',
              'specificity multiplets',
              'is cell%s' %is_cell_gex,
              'is viable cell',
              'icon','icon matching HLA']
    palette = ['grey','#ccad60','yellow','#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','black','#d6b4fc','violet'] #'#253494'
    
    #flt_to_remove = list()
    #for i, flt in enumerate(filterings):
    #    remaining_gems = sum(flt)
    #    print(labels[i], remaining_gems)
    #    if remaining_gems == 0:
    #        print('removed:', labels[i])
    #        flt_to_remove.append(i)
    
    flt_to_remove = list()
    filterings = [idx0]
    for i, flt in enumerate([idx1, idxX, idx3, idx2, idx4, idx5, idx6, idx8], start=1):
        remaining_gems = sum(filterings[-1] & flt)
        print(labels[i], remaining_gems)
        if remaining_gems > 0:
            filterings.append((filterings[-1] & flt))
        else:
            print('removed:', labels[i])
            flt_to_remove.append(i)
    # ICON HACK
    filterings.append(idx9)
    filterings.append((idx9 & idx10))
            
    for i in flt_to_remove[::-1]:
        del labels[i]
        del palette[i]
else:
    print('filterset name unknown')

##########################################################
#                       Output prep                      #
##########################################################
print('prepping output')
peptides = [peptide]
plateaus = ['intra','inter']
random_sampling = [random_sample] #np.arange(30)
clonotypes = df.ct.unique()

indexes = len(peptides)*len(plateaus)*len(labels)*len(random_sampling)*len(clonotypes)

# Initiate dataframe
test_df = pd.DataFrame(columns=['score','cdr3s'], index=np.arange(indexes))
test_df.index = pd.MultiIndex.from_product([plateaus, random_sampling, labels, peptides, clonotypes],
                                           names=['plateau', 'rnd_sample','filtering','peptide','ct'])


##########################################################
#                   Compare per peptide                  #
##########################################################
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

for data_type, flt in zip(labels, filterings):
    sdf = df[flt]
    sdf = sdf[keep_max_conc(sdf)]
    
    group = sdf[sdf.peptide_HLA == peptide].copy() #sdf.groupby('peptide_HLA', sort=False).get_group(peptide)
    #group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB'], inplace=True) # new! group hasn't be reduced earlier
    g = group[group.cdr3_TRA.isin(idx_tra[idx_tra_sub]) &
              group.cdr3_TRB.isin(idx_trb[idx_trb_sub])]
    g.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB'], inplace=True)
    
    if len(g) == 1:
        continue

    inter_chains = sdf.loc[sdf.peptide_HLA != peptide, ['cdr3_TRA', 'cdr3_TRB']] 

    # OBS! Maybe get the pairs directly instead of zipping them?
    # Make sure you get true pairs and not just random pairs?!
    #cdr3_TRAs = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRA.values
    #cdr3_TRBs = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRB.values
    
    cdr3_TRAs = g.cdr3_TRA.values #[g.cdr3_TRA.isin(idx_tra[idx_tra_sub])]
    cdr3_TRBs = g.cdr3_TRB.values #[g.cdr3_TRB.isin(idx_trb[idx_trb_sub])]
    
    print(data_type, filter_set, peptide, len(group), len(g), len(cdr3_TRAs), len(cdr3_TRBs))

    #assert len(cdr3_TRAs) == len(cdr3_TRBs) == len(group.loc[:,['cdr3_TRA','cdr3_TRB']].drop_duplicates())

    #for random_sample in random_sampling:
    intra_score_peptide = list()
    inter_score_peptide = list()
    for index, (a,b) in enumerate(zip(cdr3_TRAs, cdr3_TRBs)):
        print(index, len(g), a, b, peptide, data_type)
        intra = get_intra_similarity(cdr3_TRAs, cdr3_TRBs)
        inter = get_inter_similarity(cdr3_TRAs, cdr3_TRBs)

        intra_score_peptide.append(intra['score'])
        inter_score_peptide.append(inter['score'])

    for pl, lst in zip(['intra','inter'], [intra_score_peptide, inter_score_peptide]):
        tmp = add_data(lst, peptide, pl, data_type, random_sample)
        test_df.loc[tmp.index, ['score','cdr3s']] = tmp


test_df.reset_index(inplace=True)
test_df.score = test_df.score.astype(float)
#test_df['label'] = np.where(test_df.plateau == 'intra', 2, 0)
test_df.dropna(inplace=True) # new addition

##########################################################
#               Compute AUC for each sample              #
##########################################################
#for s, s_grp in test_df.dropna().groupby('rnd_sample'):
#    print(s)
#    for f, f_grp in s_grp.groupby('filtering'):
#        fpr, tpr, _ = roc_curve(f_grp.label, f_grp.score, pos_label=2)
#        test_df.loc[f_grp.index, 'AUC'] = auc(fpr, tpr)
#        test_df.loc[f_grp.index, 'AUC 0.1'] = auc(fpr[fpr < 0.1], tpr[fpr < 0.1])/ 0.1
#
#auc_df = test_df.dropna().melt(id_vars=['rnd_sample', 'filtering'], value_vars=['AUC','AUC 0.1'])
#auc_df['palette'] = auc_df.filtering.map(dict(zip(labels,palette)))

##########################################################
#                  Write output to file                  #
##########################################################
test_df.to_csv(SIM_OUTPUT, index=False, header=False)
#auc_df.to_csv(AUC_OUTPUT, index=False)



