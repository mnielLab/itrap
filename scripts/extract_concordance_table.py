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

# Data
INPUT = snakemake.input.valid_df #'/data/exp9.2_CAT_IONTORRENT_KMA_AKB/reports/valid_ct.csv'
OUTPUT = snakemake.output[0] #'/data/exp9.2_CAT_IONTORRENT_KMA_AKB/reports/concordance_gridsearch.relative.csv' #.relative

D_UMI_MHC = float(snakemake.params[0]) #int()
print('delta UMI MHC', D_UMI_MHC)

# Main
df = pd.read_csv(INPUT, converters=converters).fillna(value={"umi_count_mhc": 0, "delta_umi_mhc": 0,
                                                             "umi_count_mhc_rel":0,
                                                             "umi_count_cd8": 0, "delta_umi_cd8": 0,
                                                             "umi_count_TRA": 0, "delta_umi_TRA": 0,
                                                             "umi_count_TRB": 0, "delta_umi_TRB": 0})

#value_bin = df[df.valid_ct == True].ct.unique() # sign. peptides & HLA
value_bin = df[~df.ct_pep.isna()].ct.unique() # sign. peptides

# Set range of thresholds
umi_count_TRA_l = np.arange(0, df.umi_count_TRA.quantile(0.4, interpolation='higher')) #2
delta_umi_TRA_l = 2**np.linspace(-0.4,1.5,10) #np.arange(0, df.delta_umi_TRA.quantile(0.4, interpolation='higher')) #3
umi_count_TRB_l = np.arange(0, df.umi_count_TRB.quantile(0.4, interpolation='higher')) #5
delta_umi_TRB_l = 2**np.linspace(-0.4,1.5,10) #np.arange(0, df.delta_umi_TRB.quantile(0.4, interpolation='higher')) #10
umi_count_mhc_l = np.arange(1, df.umi_count_mhc.quantile(0.5, interpolation='higher')) #25
delta_umi_mhc_l = [D_UMI_MHC] #np.arange(0, df.delta_umi_mhc.quantile(0.5, interpolation='higher')) #10
umi_relat_mhc_l = [0] #[D_UMI_MHC/10000]
# umi_relat_mhc_l = np.linspace(0,0.2,20) #

#single_barcode_mhc_l = all_combinations(df.single_barcode_mhc.unique())
#exclude_specificity_singlet_l = [True, False]
#tcr_category_l = all_combinations(df.tcr_category.unique())

#print('umi_count_mhc_l', df.umi_count_mhc.quantile(0.5, interpolation='higher'))
#print('delta_umi_mhc_l', df.delta_umi_mhc.quantile(0.5, interpolation='higher'))

observations = (len(umi_count_TRA_l) *
                len(delta_umi_TRA_l) *
                len(umi_count_TRB_l) *
                len(delta_umi_TRB_l) *
                len(umi_count_mhc_l) *
                len(delta_umi_mhc_l) *
                len(umi_relat_mhc_l))
#print(observations, 'total observations')

features = ['accuracy','n_matches','n_mismatches',
            'trash_gems','ratio_retained_gems','trash_cts','ratio_retained_cts','avg_conc', 
            #'ct','gems_per_ct','ct_pep',
            'umi_count_mhc','umi_relat_mhc_l','delta_umi_mhc',
            'umi_count_TRA','delta_umi_TRA',
            'umi_count_TRB','delta_umi_TRB',
            'trash_gems_total','ratio_retained_gems_total','trash_cts_total','ratio_retained_cts_total'] #'umi_count_TRA','delta_umi_TRA',
table = pd.DataFrame(columns=features) #index=np.arange(observations + 1), 

N_total_gems = len(df)
N_total_tcrs = len(df.ct.unique())
n_total_gems = len(df[df.ct.isin(value_bin)])
n_total_tcrs = len(value_bin)

i = -1
for uca in umi_count_TRA_l:
    for dua in delta_umi_TRA_l:
        for ucb in umi_count_TRB_l:
            for dub in delta_umi_TRB_l:
                #for tc in tcr_category_l:
                for ucm in umi_count_mhc_l:
                    for urm in umi_relat_mhc_l:
                        for dum in delta_umi_mhc_l:
                            #for sbm in single_barcode_mhc_l:
                            #for ess in exclude_specificity_singlet_l:
                            i += 1
                            filter_bool = ((df.umi_count_TRA >= uca) & #| #
                                           (df.delta_umi_TRA >= dua) & #| #
                                           (df.umi_count_TRB >= ucb) & #| #
                                           (df.delta_umi_TRB >= dub) & #| #
                                           (df.umi_count_mhc >= ucm) & #| #
                                           (df.delta_umi_mhc >= dum) & #| #
                                           (df.umi_count_mhc_rel >= urm)) #(df.umi_count_mhc_rel >= urm) &
                            #(df.umi_count_TRA >= uca) & (df.delta_umi_TRA >= dua) &
                            
                            flt = df[filter_bool & df.ct.isin(value_bin)].copy()
                            flt = calc_binding_concordance(flt, 'ct')

                            n_gems = len(flt)
                            n_tcrs = len(flt.ct.unique())

                            conc = flt.binding_concordance.mean()
                            
                            assert not flt.pep_match.isna().any(), 'Make sure flt only contains value cts'
                            assert n_gems == len(flt.pep_match.dropna())
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
                                            conc, ucm, urm, dum, uca, dua, ucb, dub,
                                            G_trash, G_ratio, T_trash, T_ratio) #uca, dua, 

                            if i % 1000 == 0:
                                print(f'{round(i/observations * 100, 2)}% done')

table.to_csv(OUTPUT, index=False)
