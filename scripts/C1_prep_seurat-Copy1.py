#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import pandas as pd
from scipy.io import mmread
import os
from pathlib import Path

# ## Input data
HTO = snakemake.input.hto
BRC = snakemake.input.brc 
GEM = snakemake.input.gem
MTX = snakemake.input.mtx 
DF = snakemake.input.df #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_CAT_IONTORRENT_KMA_AKB/tables/tcr_barcode.csv' #

# ## Output data
output = snakemake.output[0] 

# ## Load
df = pd.read_csv(DF)
gems = df.dropna(subset=['ct','peptide_HLA'], how='any').gem

hto_df = pd.read_excel(HTO, sheet_name=1, names=['barcode','sample_id','HLA_A','HLA_B','HLA_C','comment']) #
htos = hto_df.barcode.astype(str)

g = np.loadtxt(GEM, dtype='U36')
b = pd.read_csv(BRC, sep='\t', header=None, names=['barcode','name','feature']) #
m = mmread(MTX)
m = pd.DataFrame.sparse.from_spmatrix(m, index=b.barcode, columns=g)

m = m.loc[htos, m.columns.isin(gems)]

if m.shape[0] == 0:
    Path(output).touch()
else:
    print(snakemake.params[0])
    Path(snakemake.params[0]).touch()
    # ## Write to rds
    import pyreadr
    pyreadr.write_rds(output, m.sparse.to_dense())


