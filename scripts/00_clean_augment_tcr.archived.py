#!/usr/bin/env python
# coding: utf-8

# # OBS! Is now a python script

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import json


# # Args
keep_only_is_cells = False
keep_only_high_confidence = False
keep_only_full_length = True
keep_only_productive = True
keep_only_unamibiguous_gems = False


# # Input
#file = "/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_TCR/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv"#"src/MS_2x150_TCR_190501/all_contig_annotations.csv" # /Volumes/tuba "src/all_contig_annotations.csv" #
INPUT = "/home/tuba/herpov/netTCR_public_data/data/donor4/all_contig_annotations.csv" #snakemake.input[0] #"/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_TCR/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv"

# get the software package from input
regex = re.search('\/exp\d\.?(2\.2)_TCR', INPUT)
if regex is not None:
	SOFTWARE = 'v2.2'
else:
	SOFTWARE = 'newer'
if 'donor' in INPUT:
    SOFTWARE = 'v3'

# # Output
#file_out = "/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_TCR/augmented/tcr.clean.augmented.csv"
OUTPUT = "/home/tuba/herpov/netTCR_public_data/data/donor4/tcr.clean.augmented.csv" #snakemake.output[0] #"/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_TCR/augmented/tcr.clean.augmented.old.csv" #
report = "/home/tuba/herpov/netTCR_public_data/data/donor4/tcr.clean.augmented.txt" #snakemake.output[1] #"/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_TCR/augmented/tcr.clean.augmented.old.txt" #

# # Load
df = pd.read_csv(INPUT) #, header=0

# ## Rename
df.rename(columns={"barcode" : "gem"}, inplace=True)
df.rename(columns={"raw_clonotype_id" : "clonotype"}, inplace=True)
df.rename(columns={"umis" : "umi_count"}, inplace=True)
df.rename(columns={"reads" : "read_count"}, inplace=True)

rprt = {'total_tcr_gems': list(set(df.gem))}

# ## Filter data
flt_report = 'filter'
if keep_only_is_cells:
    df = df[df.is_cell == True]
    flt_report += '_is_cell'
if keep_only_high_confidence:
    df = df[df.high_confidence == True]
    flt_report += '_high_confidence'
if keep_only_full_length:
    df = df[df.full_length == True]
    flt_report += '_full_length'
if keep_only_productive:
    if SOFTWARE == 'newer':
        df = df[df.productive == True]
    else:
        df = df[df.productive == 'True']
    flt_report += '_productive'
if keep_only_unamibiguous_gems:
    df = df.groupby(['gem', 'chain']).filter(lambda x: len(x) == 1)
    flt_report += '_unambiguous'

rprt[flt_report] = list(set(df.gem))

# ## Augment by chain
df['chain_count'] = df.groupby(['gem', 'chain']).contig_id.transform('size')

tra_df = df[(df.chain == "TRA")].copy()
trb_df = df[(df.chain == "TRB")].copy()

rprt['tra_gems'] = list(set(tra_df.gem))
rprt['trb_gems'] = list(set(trb_df.gem))

tra_df.sort_values(by=['gem', 'umi_count'], inplace=True)
trb_df.sort_values(by=['gem', 'umi_count'], inplace=True)

# Replacing genes annotated with None with ''
tra_df['genes'] = tra_df.replace([None], ['']).apply(lambda x: ';'.join(x[['v_gene','j_gene','c_gene']]), axis=1)
trb_df['genes'] = trb_df.replace([None], ['']).apply(lambda x: ';'.join(x[['v_gene','d_gene','j_gene','c_gene']]), axis=1)

def annotate_umi_lst(df):
    # I have to sum the UMI counts because in some GEMs the same CDR3 will be represented with different other stats
    dct = df.groupby(['gem','cdr3']).umi_count.apply(sum).to_frame().groupby('gem').umi_count.apply(lambda x: sorted(np.array(x))).to_dict()
    return df.gem.map(dct)

def annotate_cdr3_lst(df):
    # I use unique instead of listing because the same CDR3 can occur multiple times in a GEM.
    dct = df.groupby(['gem']).cdr3.unique().to_dict()
    return df.gem.map(dct)

def annotate_genes_lst(df):
    # I use unique instead of listing because the same CDR3 can occur multiple times in a GEM.
    dct = df.groupby(['gem']).genes.unique().to_dict()
    return df.gem.map(dct)

def annotate_single(df):
    return df.umi_count_lst.apply(lambda x: True if len(x)==1 else False)

def annotate_delta_umi(df):
    #def calc_delta(x):
    #    if len(x) == 1:
    #        return 100
    #    else:
    #        return int((x[-1]-x[-2])/float(x[-1])*100)
    def calc_delta(x):
        if len(x) == 1:
            return x[-1]/0.25
        elif len(x) == 0:
            return 0
        else:
            return (x[-1])/(x[-2]+0.25)
    return df.umi_count_lst.apply(calc_delta)

tra_df['umi_count_lst'] = annotate_umi_lst(tra_df)
trb_df['umi_count_lst'] = annotate_umi_lst(trb_df)

tra_df['cdr3_lst'] = annotate_cdr3_lst(tra_df)
trb_df['cdr3_lst'] = annotate_cdr3_lst(trb_df)

tra_df['genes_lst'] = annotate_genes_lst(tra_df)
trb_df['genes_lst'] = annotate_genes_lst(trb_df)

tra_df['single'] = annotate_single(tra_df)
trb_df['single'] = annotate_single(trb_df)

tra_df['delta_umi'] = annotate_delta_umi(tra_df)
trb_df['delta_umi'] = annotate_delta_umi(trb_df)

# Assert that the annotated CDR3 lst is as long as the UMI lst
assert tra_df.apply(lambda x: True if len(x.umi_count_lst) == len(x.cdr3_lst) else False, axis=1).all()


# ### Keep chain with highest UMI count and merge
tra_df.drop_duplicates(subset=['gem'], keep='last', inplace=True)
trb_df.drop_duplicates(subset=['gem'], keep='last', inplace=True)

# It is not necessary to merge on clonotype, only to keep a single column of clonotype instead of two
tcr_df = pd.merge(tra_df, trb_df, how='outer', on=['gem','clonotype'], suffixes=('_TRA', '_TRB'))

assert tcr_df.shape[0] == tcr_df.gem.unique().shape[0]
assert tcr_df.gem.unique().shape[0] == df[df.chain.isin(['TRA','TRB'])].gem.unique().shape[0], print('GEMs in new df: %d | GEMs in old df: %d' %(tcr_df.gem.unique().shape[0], df.gem.unique().shape[0]))


# ## Augment by GEM
tcr_df['num_clonotype'] = pd.to_numeric(tcr_df['clonotype'].fillna('None').str.split('clonotype').str[1],
                                        errors='coerce').replace(np.nan, 0, regex=True).astype(int)

tcr_df['single_chain_only'] = tcr_df[['chain_TRA', 'chain_TRB']].isna().any(axis=1)
tcr_df['umi_count_tcr'] = tcr_df.umi_count_TRA.fillna(0) + tcr_df.umi_count_TRB.fillna(0)
tcr_df['cdr3_comb'] = tcr_df.cdr3_TRA.fillna('') + tcr_df.cdr3_TRB.fillna('')

def define_tcr_categories(row):
    if (row['single_TRA'] == True) & (row['single_TRB'] == True):
        return 'unique chains'
    if row['single_chain_only']:
        return 'missing chain'
    else:
        return 'multiple chains'

tcr_df['tcr_category'] = tcr_df.apply(lambda row: define_tcr_categories(row), axis=1)
tcr_df['no_filtration'] = True
tcr_df['exclude_single-chain_TCRs'] = tcr_df.apply(lambda row: True if (row.tcr_category == 'unique chains') or (row.tcr_category == 'multiple chains') else False, axis=1)
tcr_df['exclude_ambiguous_and_single-chain_TCRs'] = np.where((tcr_df.tcr_category == 'unique chains'), True, False)

# ### Assign new clonotype
# I assume that a clonotype consist of a specific set of V(D)J genes from both chains and a specific combination of alpha and beta CDR3s.

def assign_clonotype():
    clonotype_variables = ['v_gene_TRA','j_gene_TRA',
                           'v_gene_TRB','j_gene_TRB',
                           'cdr3_comb']
    tcr_df.loc[:, clonotype_variables] = tcr_df.loc[:, clonotype_variables].fillna('unknown')
    new_clonotype = tcr_df.groupby(clonotype_variables).gem.unique().to_frame()
    new_clonotype['n_gems'] = new_clonotype.gem.apply(len)
    new_clonotype.sort_values(by='n_gems', ascending=False, inplace=True)
    dct = new_clonotype.to_dict()['gem']
    for i, k in enumerate(dct.keys(), start=1): 
        dct[k] = i
    return tcr_df.set_index(clonotype_variables).index.map(dct)

tcr_df['ct'] = assign_clonotype()
assert tcr_df.ct.isna().any() == False


# ## Write data
new_column_order = ['gem', 'clonotype', 'num_clonotype', 'ct', 'genes_TRA', 'genes_TRB', 'genes_lst_TRA', 'genes_lst_TRB',
                    'length_TRA', 'cdr3_TRA', 'umi_count_TRA', 'umi_count_lst_TRA', 'delta_umi_TRA', 'cdr3_lst_TRA', 'chain_count_TRA','single_TRA',
                    'length_TRB', 'cdr3_TRB', 'umi_count_TRB', 'umi_count_lst_TRB', 'delta_umi_TRB', 'cdr3_lst_TRB', 'chain_count_TRB','single_TRB',
                    'single_chain_only', 'umi_count_tcr', 'cdr3_comb', 'v_gene_TRA','j_gene_TRA', 'v_gene_TRB','j_gene_TRB', 'tcr_category', 'no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs']

tcr_df[new_column_order].to_csv(OUTPUT, index=False)

with open(report, 'w') as outfile:
    json.dump(rprt, outfile)


