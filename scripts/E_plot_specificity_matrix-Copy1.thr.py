#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import itertools
from ast import literal_eval
import re
from scipy import stats
from random import sample
#import scipy.stats as stats
import random
import re
import statistics

from D_plot_specificity_matrix_utils import (peptides_per_gem,
                                             peptide_per_clonotype_by_gem_size,
                                             multiple_peptides_per_gem_w_filtering,
                                             peptide_per_clonotype_by_umi_counts,
                                             mhc_read_count_per_clonotype,
                                             tcr_read_count_per_clonotype)

plt.style.use('ggplot')

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
              'umi_count_lst_mhc': literal_converter, #literal_eval,
              'umi_count_lst_TRA': literal_converter,'umi_count_lst_TRB': literal_converter,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'cdr3_lst_TRB': cdr3_lst_converter,
              'HLA_lst_mhc': cdr3_lst_converter,'HLA_cd8': HLA_cd8_converter} #

#pattern = re.compile('\/eval_clonotypes\/specificity_matrix\/(\w+)\.(\w+)\.(\w+)\.')
#/eval_clonotypes/specificity_matrix/{plot_type}.{data_type}.pdf

def get_filtered_data(df, opt_thr):
    # Convert series of opt_thr to an exacutable and execute
    selection = eval(' & '.join(['(df.%s >= %d)'%(k,v) for k,v in opt_thr.iterrows()]))
    assert selection.any()
    return selection

def notnan(x):
    return x == x

def get_matching_hashing(df):
    return df.apply(lambda row: row.peptide_HLA.split()[-1] in row.HLA_cd8 if (notnan(row.peptide_HLA) & notnan(row.HLA_cd8)) else False, axis=1)

def get_multiplets(df, idx1, idx2):
    tmp = df[idx1 & idx2]
    dct = tmp.groupby(['ct','peptide_HLA']).gem.count() > 1
    idx = df.set_index(['ct','peptide_HLA']).index.map(dct)
    return idx.fillna(False)

def get_complete_tcrs(df):
    return df['exclude_single-chain_TCRs']

#def get_thresholded_data(df, thrs_type, opt_thr):
#    if thrs_type == 'raw':
#        return df.any(axis=1)
#    elif thrs_type == 'thr':
#        return get_filtered_data(df, opt_thr)
#    else:
#        print('Did not receive any of the expected threshold types (raw, thr)')

def get_data(df, thrs_levl, opt_thr):
    # thrs_levl = ['1_thr_only','2_matching_hashing_only','3_specificity_multiplets_only','4_complete_tcrs_only']
    idx1 = get_filtered_data(df, opt_thr)
    idx2 = get_matching_hashing(df)
    idx3 = get_multiplets(df, idx1, idx2)
    idx4 = get_complete_tcrs(df)
    
    if thrs_levl == '1_thr_only':
        return df[idx1]
    elif thrs_levl == '2_matching_hashing_only':
        return df[idx1 & idx2]
    elif thrs_levl == '3_specificity_multiplets_only':
        return df[idx1 & idx2 & idx3]
    elif thrs_levl == '4_complete_tcrs_only':
        return df[idx1 & idx2 & idx3 & idx4]
    else:
        print('No threshold level with that name!!')

# # Input
VALID = snakemake.input.valid_df
#DATA = snakemake.input.original
THRESHOLD = snakemake.input.thresholds

# # Load
df = pd.read_csv(VALID, converters=converters) #.fillna('')
#df = pd.read_csv(DATA, converters=converters)
df.rename(columns={'rank':'epitope_rank'}, inplace=True)
df.fillna(value={"umi_count_mhc": 0, "delta_umi_mhc": 0, "umi_count_mhc_rel":0,
                 "umi_count_cd8": 0, "delta_umi_cd8": 0,
                 "umi_count_TRA": 0, "delta_umi_TRA": 0,
                 "umi_count_TRB": 0, "delta_umi_TRB": 0}, inplace=True) #df.fillna(0, inplace=True)

valid_ct = df[df.valid_ct == True].ct.unique() # sign. peptides & HLA
#valid_ct = df[~df.ct_pep.isna()].ct.unique() # sign. peptides

# # PLOT
#for filename in snakemake.output.plots:
filename = snakemake.output.plots

plot_type = snakemake.params.plot_type
thrs_levl = snakemake.params.thrs_levl

opt_thr = pd.read_csv(THRESHOLD, index_col=0, header=None)
opt_thr.dropna(inplace=True) #.fillna(0, inplace=True) # Hacks

grid = get_data(df, thrs_levl, opt_thr)

if plot_type == 'peptide_per_clonotype_by_gem_size':
    print(plot_type)
    print(filename)
    peptide_per_clonotype_by_gem_size(grid, show=False, save_tuba=filename, save_report=False)

if plot_type == 'peptides_per_gem':
    print(plot_type)
    print(filename)
    peptides_per_gem(grid, save_tuba=filename, show=False, save_report=False)

if plot_type == 'multiple_peptides_per_gem_w_filtering':
    print(plot_type)
    print(filename)
    multiple_peptides_per_gem_w_filtering(grid, show=False, save_tuba=filename, save_report=False)

if plot_type == 'mhc_read_count_per_clonotype':
    print(plot_type)
    print(filename)
    mhc_read_count_per_clonotype(grid, show=False, save_tuba=filename, save_report=False)

if plot_type == 'tcr_read_count_per_clonotype':
    print(plot_type)
    print(filename)
    tcr_read_count_per_clonotype(grid, show=False, save_tuba=filename, save_report=False)

