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

#converters = {'peptide_HLA_lst': peptide_hla_converter,
#              'umi_count_lst_mhc': literal_eval,
#              'umi_count_lst_TRA': literal_converter,'umi_count_lst_TRB': literal_converter,
#              'cdr3_lst_TRA': cdr3_lst_converter,
#              'cdr3_lst_TRB': cdr3_lst_converter,
#              'HLA_lst_mhc': cdr3_lst_converter,'HLA_cd8': HLA_cd8_converter} #

converters = {'peptide_HLA_lst': peptide_hla_converter,
              'umi_count_lst_mhc': literal_eval,
              'umi_count_lst_cd8': literal_converter,
              'umi_count_lst_TRA': literal_converter,'umi_count_lst_TRB': literal_converter,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'cdr3_lst_TRB': cdr3_lst_converter,
              'HLA_lst_mhc': cdr3_lst_converter,
              'HLA_pool_cd8':cdr3_lst_converter,
              'HLA_cd8': HLA_cd8_converter,
              'HLA_lst_cd8':literal_converter,'sample_id_lst':literal_converter} #

def plot_grid(grid, opt_thr_idx, output):
    x_min = grid.index.max()*0.01
    x_max = grid.index.max()

    fig = plt.figure(figsize=(16,9))
    ax = plt.gca()
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)

    x = grid.index
    y = grid.accuracy
    plt.scatter(x,y, label='Accuracy', marker='.')

    x = grid.index
    y = grid.ratio_retained_gems
    plt.scatter(x,y, label='Retained GEMs', marker='.')

    x = grid.index
    y = grid.mix_mean
    plt.scatter(x,y, label='Mix mean', marker='.')

    acc, rat = grid.loc[opt_thr_idx, ['accuracy','ratio_retained_gems']]

    plt.hlines(y=[acc, rat], xmin=-x_min, xmax=opt_thr_idx, colors='grey', linestyles='--')
    plt.vlines(x=opt_thr_idx, ymin=0.05, ymax=acc, colors='grey', linestyles='--')

    t = ', '.join(grid.loc[opt_thr_idx, ['umi_count_mhc','delta_umi_mhc',
                                         'umi_count_TRA','delta_umi_TRA',
                                         'umi_count_TRB','delta_umi_TRB']].astype(int).astype(str).to_list())
    plt.text(opt_thr_idx, 0.04, t, ha='center', va='top')
    plt.text(-0.01, acc, str(acc), ha='right', va='center', transform=trans)
    plt.text(-0.01, rat, str(rat), ha='right', va='center', transform=trans)

    plt.xlim(-x_min, x_max + x_min)
    plt.ylim(-0.01, 1.01)
    plt.legend(bbox_to_anchor=(0.99, 0.5), loc='center right', frameon=False)
    plt.xlabel('Grid index')
    #plt.ylabel('Fraction')

    for f in output:
        plt.savefig(f, bbox_inches='tight')

# # Input
VALID = snakemake.input.valid
GRIDS = snakemake.input.grids

# # Output
OUTPUT = snakemake.output.plots
THRESHOLD = snakemake.output.table

# # Load
valid_df = pd.read_csv(VALID, converters=converters).fillna('')

tmp = list()
for filename in GRIDS:
    tmp.append(pd.read_csv(filename))
grid = pd.concat(tmp)

# # Main
n = 2
grid['mix_mean'] = (grid.accuracy * n + grid.ratio_retained_gems)/(n+1)

# Set index according to sorting so that plot will look nicer
grid.sort_values(by=['accuracy','ratio_retained_gems'], inplace=True)
grid.reset_index(drop=True, inplace=True)

optimal_thresholds = (grid
                      .sort_values(by=['mix_mean', 'accuracy','ratio_retained_gems', #'umi_count_mhc_rel',
                                       'umi_count_mhc', 'delta_umi_mhc','umi_count_TRB','delta_umi_TRB'],
                                   ascending=[True, True, True, False, False, False, False])
                      .tail(20))

# Get index of optimal threshold
opt_thr_idx = optimal_thresholds.tail(1).index[0]
# Get threshold values
opt_thr = optimal_thresholds.loc[opt_thr_idx, ['umi_count_mhc','delta_umi_mhc', 'umi_count_mhc_rel',
                                               'umi_count_cd8','delta_umi_cd8',
                                               'umi_count_TRA','delta_umi_TRA',
                                               'umi_count_TRB','delta_umi_TRB']]

print(opt_thr_idx, opt_thr)
opt_thr.to_csv(THRESHOLD)
plot_grid(grid, opt_thr_idx, OUTPUT)

# Combine thr to string for output-names: '14.2.0.5.0.0'
# thr_name = '.'.join(opt_thr.astype(int).astype(str).values)
