#!/bin/bash/env python

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from ast import literal_eval

plt.style.use('ggplot')

import os
import sys
from D_plot_specificity_matrix_utils import (peptides_per_gem,
                                             peptide_per_clonotype_by_gem_size,
                                             multiple_peptides_per_gem_w_filtering,
                                             peptide_per_clonotype_by_umi_counts,
                                             mhc_read_count_per_clonotype,
                                             tcr_read_count_per_clonotype)

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

converters = {'peptide_HLA_lst': peptide_hla_converter,
              'epitope_lst': epitope_converter,
              'umi_count_lst_mhc': literal_eval,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'HLA_lst_mhc': cdr3_lst_converter,
              'HLA_cd8': HLA_cd8_converter}

def get_filtered_data(df, opt_thr):
    # Convert series of opt_thr to an exacutable and execute
    selection = eval(' & '.join(['(df.%s >= %d)'%(k,v) for k,v in opt_thr.iterrows()]))
    assert selection.any()
    return selection

def get_thresholded_data(df, thrs_type, opt_thr):
    if thrs_type == 'raw':
        return df.any(axis=1)
    elif thrs_type == 'thr':
        return get_filtered_data(df, opt_thr)
    else:
        print('Did not receive any of the expected threshold types (raw, thr)')

def get_data(df, data_type, thrs_type, opt_thr):
    filtering = get_thresholded_data(df, thrs_type, opt_thr)
    if data_type == 'train':
        return df[df.ct.isin(valid_ct) & filtering].copy()
    elif data_type == 'test':
        return df[~df.ct.isin(valid_ct) & filtering].copy()
    elif data_type == 'total':
        return df[filtering].copy()
    else:
        print('Did not receive any of the expected data types (train, test, total)')

##################################################################################################
#                                              Input                                             #
##################################################################################################
VALID = snakemake.input.valid_df
THRESHOLD = snakemake.input.thresholds

clonotype_fmt = 'ct'
umi_delta = 0
level = 'no_filtration'
plot_type = snakemake.params.plot_type
data_type = snakemake.params.data_type
thrs_type = snakemake.params.thrs_type

print(plot_type, data_type, thrs_type)

##################################################################################################
#                                              Load                                              #
##################################################################################################
df = pd.read_csv(VALID, converters=converters, low_memory=False)
df.rename(columns={'rank':'epitope_rank'}, inplace=True)
df.fillna(value={"umi_count_mhc": 0, "delta_umi_mhc": 0, "umi_count_mhc_rel":0,
                 "umi_count_cd8": 0, "delta_umi_cd8": 0,
                 "umi_count_TRA": 0, "delta_umi_TRA": 0,
                 "umi_count_TRB": 0, "delta_umi_TRB": 0}, inplace=True)

opt_thr = pd.read_csv(THRESHOLD, index_col=0, header=None)
opt_thr.dropna(inplace=True) #.fillna(0, inplace=True) # Hacks

valid_ct = df[df.valid_ct == True].ct.unique() # sign. peptides & HLA
#valid_ct = df[~df.ct_pep.isna()].ct.unique() # sign. peptides

grid = get_data(df, data_type, thrs_type, opt_thr)
grid.dropna(subset=['sample_id'], inplace=True)

output_path = os.path.dirname(os.path.abspath(snakemake.output[0]))

##################################################################################################
#                                              Plot                                              #
##################################################################################################
for sample, grp in grid.groupby('sample_id'):
    if grp.empty:
        continue
    for extension in ['.pdf', '.png']:
        filename = output_path + '/' + str(int(sample)) + extension
        print(filename)

        if plot_type == 'peptide_per_clonotype_by_gem_size':
            print(plot_type)
            peptide_per_clonotype_by_gem_size(grp, show=False, save_tuba=filename, save_report=False)

        if plot_type == 'peptides_per_gem':
            print(plot_type)
            peptides_per_gem(grid, save_tuba=filename, show=False, save_report=False)

        if plot_type == 'multiple_peptides_per_gem':
            print(plot_type)
            print('Not implemented')

        if plot_type == 'multiple_peptides_per_gem_w_filtering':
            print(plot_type)
            multiple_peptides_per_gem_w_filtering(grp, show=False, save_tuba=filename, save_report=False)

        if plot_type == 'mhc_read_count_per_clonotype':
            print(plot_type)
            mhc_read_count_per_clonotype(grp, show=False, save_tuba=filename, save_report=False)

        if plot_type == 'tcr_read_count_per_clonotype':
            print(plot_type)
            tcr_read_count_per_clonotype(grp, show=False, save_tuba=filename, save_report=False)

