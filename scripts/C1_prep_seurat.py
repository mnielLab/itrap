#!/usr/bin/env python
# coding: utf-8

# In[1]:

#import pyreadr
import numpy as np
import pandas as pd
from scipy.io import mmread


# ## Input data
HTO = snakemake.input.hto
BRC = snakemake.input.brc 
GEM = snakemake.input.gem
MTX = snakemake.input.mtx 
DF = snakemake.input.df

# ## Output data
output = snakemake.output[0] #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_CAT_IONTORRENT_KMA_AKB/tables/tcr_barcode.csv' #

# ## Load
df = pd.read_csv(DF)
gems = df.dropna(subset=['ct','peptide_HLA'], how='any').gem

hto_df = pd.read_excel(HTO, sheet_name=1, names=['barcode','sample_id','HLA_A','HLA_B','HLA_C','comment'])
htos = hto_df.barcode

g = np.loadtxt(G, dtype='U36')
b = pd.read_csv(B, sep='\t', header=None, names=['barcode','name','feature'])
m = mmread(M)
matrix = pd.DataFrame(m.toarray(), index=b.barcode, columns=g)
    
    


# ## Write to rds
#pyreadr.write_rds(output, matrix.loc[htos, gems])


