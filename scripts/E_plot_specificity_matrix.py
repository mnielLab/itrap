#!/bin/bash/env python
#https://github.com/ipython/ipynb
#https://stackoverflow.com/questions/34478398/import-local-function-from-a-module-housed-in-another-directory-with-relative-im

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from ast import literal_eval

plt.style.use('ggplot')

import os
import sys
from D_plot_specificity_matrix_utils import (epitope_sorter_index,
											 filtering,
                       get_nonsinglet_idxs,
											 peptides_per_gem,
											 peptide_per_clonotype_by_gem_size,
                       multiple_peptides_per_gem_w_filtering,
                                             peptide_per_clonotype_by_umi_counts,
                                             mhc_read_count_per_clonotype,
                                             tcr_read_count_per_clonotype,
                                             multiple_peptides_per_gem)
#module_path = os.path.abspath(os.path.join('../tcr-pmhc-sc-project/notebooks/'))
#if module_path not in sys.path:
#	sys.path.append(module_path)
#from ipynb.fs.full.D_plot_specificity_matrix import calc_binding_concordance, peptides_per_gem, peptide_per_clonotype_by_gem_size, peptide_per_clonotype_read_counts, mhc_read_count_per_clonotype, mhc_read_count_per_clonotype_response, mhc_read_count_per_clonotype_peptide_assayed

#def cdr3_lst_converter(x):
#    #define format of datetime
#    return x.replace("[","").replace("]","").replace("'","").split(" ")
#
#def epitope_converter(x):
#    #define format of datetime
#    return x.replace("[","").replace("]","").replace("'","").split(" ")
#
#def peptide_hla_converter(x):
#    return re.findall("\w+\s{1}\w{1}\d+", x.replace("[","").replace("]","").replace("'",""))
#
#converters = {'peptide_HLA_lst': peptide_hla_converter,
#			  'epitope_lst': epitope_converter,
#			  'umi_count_lst_mhc': literal_eval,
#			  'cdr3_lst_TRA': cdr3_lst_converter}

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

# Get first file and extract information, do filtering and ranking. Then loop through and plot:
#filename = snakemake.output.plots[0]
#print(filename)
#pattern = re.compile('\/specificity_matrix\/(\w+)\/(\w+)\/umi_delta(\d+)\/(\w+\-?\w*)\/b(\d+).t(\d+).ecs_(\w+).ess_(\w+).')
#regex = pattern.search(filename)
plot_types = snakemake.params.plot_type #regex.group(1)
clonotype_fmt = 'ct' #regex.group(2)
umi_delta = 0 #int(regex.group(3))
level = 'no_filtration' #regex.group(4)
bc_threshold = int(snakemake.params.min_brc[1:]) #int(regex.group(5))
tcr_threshold = int(snakemake.params.min_tcr[1:]) #int(regex.group(6))
ecs = literal_eval(snakemake.params.ecs[4:]) #literal_eval(regex.group(7))
ess = literal_eval(snakemake.params.ess[4:]) #literal_eval(regex.group(8))

print('BRC threshold:', bc_threshold)
print('TCR threshold:', tcr_threshold)

# Main
df = pd.read_csv(snakemake.input[0], converters=converters, low_memory=False)
original_length = len(df)

#df['epitope_rank'] = epitope_sorter_index(df)
df.rename(columns={'rank':'epitope_rank'}, inplace=True)

df = filtering(df, bc_threshold, tcr_threshold,
         exclude_clonotype_singlets=ecs, exclude_specificity_singlets=ess,
         clonotype_fmt=clonotype_fmt, filtration=level, umi_delta=umi_delta)

if (level == 'no_filtration') & (umi_delta == 0) & (bc_threshold == 0) & (tcr_threshold == 0):
    assert len(df) == original_length

# Extract args from filename
for plot_type, filename in zip(plot_types, snakemake.output.plots):
    print(plot_type)
    print(filename)
    #regex = pattern.search(filename)
    #plot_type = regex.group(1)
    #assert clonotype_fmt == regex.group(2)
    #assert umi_delta == int(regex.group(3))
    #assert level == regex.group(4)
    #assert bc_threshold == int(regex.group(5))
    #assert tcr_threshold == int(regex.group(6))
    #assert ecs == literal_eval(regex.group(7))
    #assert ess == literal_eval(regex.group(8))

    if plot_type == 'peptide_per_clonotype_by_gem_size':
        print(plot_type)
        peptide_per_clonotype_by_gem_size(df,
        								  clonotype_fmt=clonotype_fmt,
        								  bc_threshold=bc_threshold,
        								  tcr_threshold=tcr_threshold,
        								  exclude_clonotype_singlets=ecs,
        								  exclude_specificity_singlets=ess,
        								  show=False,
        								  save_tuba=filename,
        								  save_sund=False,
        								  save_report=False)

    if plot_type == 'peptides_per_gem':
        print(plot_type)
        peptides_per_gem(df,
                         clonotype_fmt=clonotype_fmt,
                         bc_threshold=bc_threshold,
                         tcr_threshold=tcr_threshold,
                         exclude_clonotype_singlets=ecs,
                         exclude_specificity_singlets=ess,
                         show=False,
                         save_tuba=filename,
                         save_sund=False,
                         save_report=False)

    if plot_type == 'multiple_peptides_per_gem':
        print(plot_type)
        multiple_peptides_per_gem(df,
                                  clonotype_fmt=clonotype_fmt,
                                  bc_threshold=bc_threshold,
                                  tcr_threshold=tcr_threshold,
                                  exclude_clonotype_singlets=ecs,
                                  exclude_specificity_singlets=ess,
                                  show=False,
                                  save_tuba=filename,
                                  save_sund=False,
                                  save_report=False)

    if plot_type == 'multiple_peptides_per_gem_w_filtering':
        print(plot_type)
        multiple_peptides_per_gem_w_filtering(df,
                                  clonotype_fmt=clonotype_fmt,
                                  bc_threshold=bc_threshold,
                                  tcr_threshold=tcr_threshold,
                                  exclude_clonotype_singlets=ecs,
                                  exclude_specificity_singlets=ess,
                                  show=False,
                                  save_tuba=filename,
                                  save_sund=False,
                                  save_report=False)

    if plot_type == 'mhc_read_count_per_clonotype':
        print(plot_type)
        mhc_read_count_per_clonotype(df,
                                  clonotype_fmt=clonotype_fmt,
                                  bc_threshold=bc_threshold,
                                  tcr_threshold=tcr_threshold,
                                  exclude_clonotype_singlets=ecs,
                                  exclude_specificity_singlets=ess,
                                  show=False,
                                  save_tuba=filename,
                                  save_sund=False,
                                  save_report=False)

    if plot_type == 'tcr_read_count_per_clonotype':
        print(plot_type)
        tcr_read_count_per_clonotype(df,
                                  clonotype_fmt=clonotype_fmt,
                                  bc_threshold=bc_threshold,
                                  tcr_threshold=tcr_threshold,
                                  exclude_clonotype_singlets=ecs,
                                  exclude_specificity_singlets=ess,
                                  show=False,
                                  save_tuba=filename,
                                  save_sund=False,
                                  save_report=False)

#pattern = re.compile('\/specificity_matrix\/(\w+)\/umi_delta(\d+)\/(\w+\-?\w*)\/(\w+)\/b(\d+).t(\d+).ecs_(\w+).ess_(\w+).')
#
#for filename in snakemake.output.report:
#    # Report lists of unique TCRs and GEMs
#    regex = pattern.search(filename)
#    clonotype_fmt = regex.group(1)
#    data_type = regex.group(4)
#    
#    if data_type == 'unique_tcrs':
#        np.savetxt(filename, df[clonotype_fmt].unique(), fmt='%s')
#    if data_type == 'unique_gems':
#        np.savetxt(filename, df.gem.unique(), fmt='%s')
