#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from scipy.io import mmread
import os
import argparse

def get_argparser():
    """Return an argparse argument parser."""
    parser = argparse.ArgumentParser(prog = 'Prepare Seurat Analysis',
                                     description = 'Convert data to and RDS for Seurat HTO analysis')
    add_arguments(parser)
    return parser

def add_arguments(parser):
    parser.add_argument('--features', required=True, help='Cellranger multi output: features.tsv.gz')
    parser.add_argument('--barcodes', required=True, help='Cellranger multi output: barcodes.tsv.gz')
    parser.add_argument('--matrix', required=True, help='Cellranger multi output: matrix.tsv.gz')
    parser.add_argument('--multimers', required=True, help='Filepath for multimer barcode specifications')
    parser.add_argument('--data', required=True, help='Filepath of merged TCR and barcode data')
    parser.add_argument('--output', required=True, help='Filepath of output data')
    
try:
    HTO = snakemake.input.hto
    BRC = snakemake.input.brc 
    GEM = snakemake.input.gem
    MTX = snakemake.input.mtx 
    DF = snakemake.input.df
    output = snakemake.output[0]
except:
    parser = get_argparser()
    args = parser.parse_args()
    
    HTO = args.multimers
    BRC = args.features
    GEM = args.barcodes
    MTX = args.matrix
    DF = args.data
    output = args.output

# ## Load
df = pd.read_csv(DF)
gems = df.gem # .dropna(subset=['ct','peptide_HLA'], how='any')

hto_df = pd.read_excel(HTO, sheet_name=1, names=['barcode','sample_id','HLA_A','HLA_B','HLA_C','comment'])
htos = hto_df.barcode.astype(str)

g = np.loadtxt(GEM, dtype='U36')
b = pd.read_csv(BRC, sep='\t', header=None, names=['barcode','name','feature']) #
m = mmread(MTX)
m = pd.DataFrame.sparse.from_spmatrix(m, index=b.barcode, columns=g)

m = m.loc[htos, m.columns.isin(gems)]

# ## Write to rds
import pyreadr
pyreadr.write_rds(output, m.sparse.to_dense())


