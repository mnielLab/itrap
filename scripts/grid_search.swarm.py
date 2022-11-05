#!/usr/bin/env python
#!/bin/bash/env python
#https://github.com/ipython/ipynb
#https://stackoverflow.com/questions/34478398/import-local-function-from-a-module-housed-in-another-directory-with-relative-im

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from ast import literal_eval
from itertools import chain, combinations

plt.style.use('ggplot')

import os
import sys
from D_plot_specificity_matrix_utils import (calc_binding_concordance)

def all_combinations(ss):
    #https://stackoverflow.com/questions/464864/how-to-get-all-possible-combinations-of-a-list-s-elements
    # Get all combinations of a list's elements
    return [list(subset) for subset in chain(*map(lambda x: combinations(ss, x), range(1, len(ss)+1)))]


def HLA_cd8_converter(x):
    # "['A0201', 'A2501', 'B0702', 'B3501', 'C0401', 'C0702']"
    return x.replace("[","").replace("]","").replace(",", "").replace("'","").split(" ")

def cdr3_lst_converter(x):
    # "['A0201' 'A0101' 'B0801']"
    return x.replace("[","").replace("]","").replace("'","").split(" ")

def epitope_converter(x):
    # "['06_1_1' '45_1_49' 'V15_A1  CMV  pp50  VTE' '45_1_3' '47_1_78'\n 'V17_B8 EBV BZLF1 (C9)']"
    return [y for y in x.replace("[","").replace("]","").replace("\n","").split("'") if (y != '') & (y != ' ')]

def peptide_hla_converter(x):
    # "['ALPGVPPV A0201' 'RQAYLTNQY A0101' 'VTEHDTLLY A0101' 'EERQAYLTNQY A0101'\n 'AMLIRDRL B0801' 'RAKFKQLL B0801' 'p1.a1 *A0101']"
    # Added: |p1.a1 \*\w\d{4}
    return re.findall("\w+\s{1}\w{1}\d+|p1.a1 p\*\w\d{4}", x.replace("[","").replace("]","").replace("\n","").replace("'",""))

def literal_converter(val):
    # replace NaN with '' and perform literal eval on the rest
    return [] if val == '' else literal_eval(val)

converters = {'peptide_HLA_lst': peptide_hla_converter,
              'umi_count_lst_mhc': literal_eval,
              'umi_count_lst_TRA': literal_converter,'umi_count_lst_TRB': literal_converter,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'cdr3_lst_TRB': cdr3_lst_converter,
              'HLA_lst_mhc': cdr3_lst_converter,'HLA_cd8': HLA_cd8_converter}

def f(a,b,c,d,e,f,g,h):
    output = list()
    for j in range(n_particles):
        flt = df[(df.umi_count_mhc >= a[j]) & (df.delta_umi_mhc >= b[j]) &
                 (df.umi_count_cd8 >= c[j]) & (df.delta_umi_cd8 >= d[j]) &
                 (df.umi_count_TRA >= e[j]) & (df.delta_umi_TRA >= f[j]) &
                 (df.umi_count_TRB >= g[j]) & (df.delta_umi_TRB >= h[j]) &
                 (df.valid_ct == True)]
        if len(flt) == 0:
            output.append(0)
        else:
            acc = flt.train_label.sum() / len(flt)
            rat = len(flt) / len(df[df.valid_ct == True])
            n = 2

            output.append(-(acc * n + rat) / (n+1))
        
    return np.array(output)

def update():
    "Function to do one iteration of particle swarm optimization"
    global V, X, pbest, pbest_obj, gbest, gbest_obj
    # Update params
    r1, r2 = np.random.rand(2)
    V = w * V + c1*r1*(pbest - X) + c2*r2*(gbest.reshape(-1,1)-X)
    X = np.round(X + V)
    obj = f(X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7])
    pbest[:, (pbest_obj >= obj)] = X[:, (pbest_obj >= obj)]
    pbest_obj = np.array([pbest_obj, obj]).min(axis=0)
    gbest = pbest[:, pbest_obj.argmin()]
    gbest_obj = pbest_obj.min()
    
    pbest_std = np.std(pbest, axis=1).mean()
    #print(gbest)
    #print(gbest_obj)
    #print(pbest_std)
    
    return pbest_std

def write_to_table(gbest):
    global table
    print(gbest)
    ucm, dum, uch, duh, uca, dua, ucb, dub = [0 if gb < 0 else gb for gb in gbest]
    urm = 0 # Hack. Not testing for relative MHC UMI counts in this script.
    print(ucm, dum, uch, duh, uca, dua, ucb, dub)
    
    filter_bool = ((df.umi_count_mhc >= ucm) & (df.delta_umi_mhc >= dum) &
                   (df.umi_count_cd8 >= uch) & (df.delta_umi_cd8 >= duh) &
                   (df.umi_count_TRA >= uca) & (df.delta_umi_TRA >= dua) &
                   (df.umi_count_TRB >= ucb) & (df.delta_umi_TRB >= dub))
    
    flt = df[filter_bool & df.ct.isin(value_bin)].copy()
    flt = calc_binding_concordance(flt, 'ct')

    n_gems = len(flt)
    n_tcrs = len(flt.ct.unique())

    conc = flt.binding_concordance.mean()

    n_mat = flt.pep_match.sum() # sign. peptides only
    #n_mat = flt.train_label.sum() # sign. peptides & HLA
    n_mis = n_gems - n_mat

    g_trash = n_total_gems - n_gems
    t_trash = n_total_tcrs - len(flt.ct.unique())
    G_trash = N_total_gems - n_gems
    T_trash = N_total_tcrs - n_tcrs

    g_ratio = round(n_gems / n_total_gems, 3)
    t_ratio = round(n_tcrs / n_total_tcrs, 3)
    G_ratio = round(n_gems / N_total_gems, 3)
    T_ratio = round(n_tcrs / N_total_tcrs, 3)

    acc = round(n_mat/n_gems, 3)

    table.loc[i] = (acc, n_mat, n_mis,
                    g_trash, g_ratio, t_trash, t_ratio,
                    conc, ucm, urm, dum, uch, duh, uca, dua, ucb, dub,
                    G_trash, G_ratio, T_trash, T_ratio) #uca, dua, 

# Data
INPUT = snakemake.input[0] #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9.2_CAT_IONTORRENT_KMA_AKB/reports/valid_ct.csv'
OUTPUT = snakemake.output[0] #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9.2_CAT_IONTORRENT_KMA_AKB/reports/concordance_gridsearch.relative.csv' #.relative

random_seed = int(snakemake.params[0])

# Main
df = pd.read_csv(INPUT, converters=converters)
df.fillna(value={"umi_count_mhc": 0, "delta_umi_mhc": 0, "umi_count_mhc_rel":0,
                 "umi_count_cd8": 0, "delta_umi_cd8": 0,
                 "umi_count_TRA": 0, "delta_umi_TRA": 0,
                 "umi_count_TRB": 0, "delta_umi_TRB": 0}, inplace=True)

value_bin = df[df.valid_ct == True].ct.unique() # sign. peptides & HLA
#value_bin = df[~df.ct_pep.isna()].ct.unique() # sign. peptides

# Set range of thresholds
umi_count_TRA_l = df.umi_count_TRA.quantile(0.5, interpolation='higher') #2
delta_umi_TRA_l = df.delta_umi_TRA.quantile(0.5, interpolation='higher') #3
umi_count_TRB_l = df.umi_count_TRB.quantile(0.5, interpolation='higher') #5
delta_umi_TRB_l = df.delta_umi_TRB.quantile(0.5, interpolation='higher') #10
umi_count_cd8_l = df.umi_count_cd8.quantile(0.5, interpolation='higher') #2
delta_umi_cd8_l = df.delta_umi_cd8.quantile(0.5, interpolation='higher') #3
umi_count_mhc_l = df.umi_count_mhc.quantile(0.5, interpolation='higher') #25
delta_umi_mhc_l = df.delta_umi_mhc.quantile(0.5, interpolation='higher') #10
#umi_relat_mhc_l = [0] #[D_UMI_MHC/10000]
# umi_relat_mhc_l = np.linspace(0,0.2,20) #

#single_barcode_mhc_l = all_combinations(df.single_barcode_mhc.unique())
#exclude_specificity_singlet_l = [True, False]
#tcr_category_l = all_combinations(df.tcr_category.unique())

#https://machinelearningmastery.com/a-gentle-introduction-to-particle-swarm-optimization/
# Hyper-parameter of the algorithm
c1 = c2 = 0.1
w = 0.8
 
# Create particles
n_particles = 1000
n_parameters = 8
np.random.seed(random_seed)

init = np.array([[umi_count_mhc_l]*n_particles,
                 [delta_umi_mhc_l]*n_particles,
                 [umi_count_cd8_l]*n_particles,
                 [delta_umi_cd8_l]*n_particles,
                 [umi_count_TRA_l]*n_particles,
                 [delta_umi_TRA_l]*n_particles,
                 [umi_count_TRB_l]*n_particles,
                 [delta_umi_TRB_l]*n_particles])
X = np.round(np.random.rand(n_parameters, n_particles) * init)
V = np.random.randn(n_parameters, n_particles) * 0.8
 
# Initialize data
pbest = X
pbest_obj = f(X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7])
gbest = pbest[:, pbest_obj.argmin()]
gbest_obj = pbest_obj.min()

features = ['accuracy','n_matches','n_mismatches',
            'trash_gems','ratio_retained_gems','trash_cts','ratio_retained_cts','avg_conc', 
            #'ct','gems_per_ct','ct_pep',
            'umi_count_mhc','umi_relat_mhc_l','delta_umi_mhc',
            'umi_count_cd8','delta_umi_cd8',
            'umi_count_TRA','delta_umi_TRA',
            'umi_count_TRB','delta_umi_TRB',
            'trash_gems_total','ratio_retained_gems_total','trash_cts_total','ratio_retained_cts_total'] #'umi_count_TRA','delta_umi_TRA',
table = pd.DataFrame(columns=features) #index=np.arange(observations + 1), 

N_total_gems = len(df) #len(original_df)
N_total_tcrs = len(df.ct.unique()) #len(original_df.ct.unique())
n_total_gems = len(df[df.ct.isin(value_bin)])
n_total_tcrs = len(value_bin)

prev_std = None
list_std = list
for i in range(100):
    curr_std = update()
    write_to_table(gbest)
    
    if curr_std == prev_std:
        list_std.append(curr_std)
    else:
        list_std = list()
        
    if len(list_std) == 10:
        break
        
    prev_std = curr_std
    
table.to_csv(OUTPUT, index=False)

#i = -1
#for uca in umi_count_TRA_l:
#    for dua in delta_umi_TRA_l:
#        for ucb in umi_count_TRB_l:
#            for dub in delta_umi_TRB_l:
#                #for tc in tcr_category_l:
#                for ucm in umi_count_mhc_l:
#                    for urm in umi_relat_mhc_l:
#                        for dum in delta_umi_mhc_l:
#                            #for sbm in single_barcode_mhc_l:
#                            #for ess in exclude_specificity_singlet_l:
#                            i += 1
#                            filter_bool = ((df.umi_count_TRA >= uca) &
#                                           (df.delta_umi_TRA >= dua) &
#                                           (df.umi_count_TRB >= ucb) &
#                                           (df.delta_umi_TRB >= dub) &
#                                           (df.umi_count_mhc >= ucm) &
#                                           (df.delta_umi_mhc >= dum) &
#                                           (df.umi_count_mhc_rel >= urm)) #(df.umi_count_mhc_rel >= urm) &
#                            #(df.umi_count_TRA >= uca) & (df.delta_umi_TRA >= dua) &
#                            
#                            flt = df[filter_bool & df.ct.isin(value_bin)].copy()
#                            flt = calc_binding_concordance(flt, 'ct')
#
#                            n_gems = len(flt)
#                            n_tcrs = len(flt.ct.unique())
#
#                            conc = flt.binding_concordance.mean()
#
#                            n_mat = flt.pep_match.sum() # sign. peptides only
#                            #n_mat = flt.train_label.sum() # sign. peptides & HLA
#                            n_mis = n_gems - n_mat
#
#                            g_trash = n_total_gems - n_gems
#                            t_trash = n_total_tcrs - len(flt.ct.unique())
#                            G_trash = N_total_gems - n_gems
#                            T_trash = N_total_tcrs - n_tcrs
#
#                            g_ratio = round(n_gems / n_total_gems, 3)
#                            t_ratio = round(n_tcrs / n_total_tcrs, 3)
#                            G_ratio = round(n_gems / N_total_gems, 3)
#                            T_ratio = round(n_tcrs / N_total_tcrs, 3)
#
#                            acc = round(n_mat/n_gems, 3)
#
#                            table.loc[i] = (acc, n_mat, n_mis, g_trash, g_ratio, t_trash, t_ratio, conc, ucm, urm, dum, uca, dua, ucb, dub, G_trash, G_ratio, T_trash, T_ratio) #uca, dua, 
#
#                            if i % 1000 == 0:
#                                print(f'{round(i/observations * 100, 2)}% done')
#
#table.to_csv(OUTPUT, index=False)
