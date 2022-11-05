#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')

# ## Input data
TCR = snakemake.input[0] #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_TCR/augmented/tcr.clean.augmented.csv' #
BRC = snakemake.input[1] #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_MHC_IONTORRENT/mapping/KMA-1t1/output/mapping.clean.AKB.augmented.gz' #

# ## Output data
output = snakemake.output[0] #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_CAT_IONTORRENT_KMA_AKB/tables/tcr_barcode.csv' #

# ## Load
tcr_df = pd.read_csv(TCR)
brc_df = pd.read_csv(BRC)
#, usecols=['query_id', 'template_id', 'gem', 'bit_score', 'alignment_length', 'tso', 'b_primer', 'anneal', 'a_primer', 'match'], sep=" ", names=["read_id", "gem", "tso", "b_primer", "anneal", "cd8_primer", "mhc_primer"]

# QC
print("Augmented barcodes")
print("GEMs: %i" %brc_df.gem.unique().shape[0])

print("Cleaned TCRs")
print("GEMs: %i" %tcr_df.gem.unique().shape[0])

# ## Merge clonotypes and barcodes
clonotype_barcode_specificity_df = pd.merge(tcr_df,
                                            brc_df,
                                            how='right', on='gem')

columns = clonotype_barcode_specificity_df.columns.to_list()
columns.remove('umi_count_lst_TRA')
columns.remove('umi_count_lst_TRB')
columns.remove('umi_count_lst_mhc')
columns.remove('umi_count_lst_cd8')
#columns.remove('read_count_lst_mhc')
#columns.remove('read_count_lst_cd8')
#columns.remove('template_lst_mhc')
#columns.remove('template_lst_cd8')
columns.remove('cdr3_lst_TRA')
columns.remove('cdr3_lst_TRB')
columns.remove('cdr3_TRA')
columns.remove('cdr3_TRB')
columns.remove('epitope_lst')
columns.remove('genes_lst_TRA')
columns.remove('genes_lst_TRB')

clonotype_barcode_specificity_df.drop_duplicates(subset=columns, inplace=True)
clonotype_barcode_specificity_df.fillna({'umi_count_mhc':0, 'delta_umi_mhc':0,
                                         'umi_count_TRA':0, 'delta_umi_TRA':0,
                                         'umi_count_TRB':0, 'delta_umi_TRB':0}, inplace=True)


# QC
print("Merged barcodes and TCRs")
print("Rows: %i" %clonotype_barcode_specificity_df.shape[0])
print("Gems: %i" %clonotype_barcode_specificity_df.gem.unique().shape[0])


# ## Write to excel
clonotype_barcode_specificity_df.to_csv(output, index=False)


