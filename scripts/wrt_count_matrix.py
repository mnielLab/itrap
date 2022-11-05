#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import os

INPUT = snakemake.input[0] #'/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run1/mhc/mapping/split/umi.tsv' #
MTX = snakemake.output.mtx
BRC = snakemake.output.brc
GEM = snakemake.output.gem

mtx, ext = os.path.splitext(MTX)

df = pd.read_csv(INPUT, header=None, names=['gem','barcode','umi','label'])

m = df.pivot(index='barcode', columns='gem', values='umi')
m.fillna(0, inplace=True)
s = csr_matrix(m)
s.eliminate_zeros()

mmwrite(mtx, s, field='integer')

pd.DataFrame([m.index, m.index, m.index]).T.to_csv(BRC, index=False, header=False, sep='\t')
pd.DataFrame(m.columns).to_csv(GEM, index=False, header=False)