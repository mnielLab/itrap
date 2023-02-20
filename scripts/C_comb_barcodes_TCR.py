#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import argparse

def get_argparser():
    """Return an argparse argument parser."""
    parser = argparse.ArgumentParser(prog = 'Merge',
                                     description = 'Merge TCR and multimer barcode data')
    add_arguments(parser)
    return parser

def add_arguments(parser):
    parser.add_argument('--tcr', required=True, help='Filepath for TCR data')
    parser.add_argument('--barcodes', required=True, help='Filepath for multimer barocde data')
    parser.add_argument('--output', required=True, help='Filepath of output data')

try:
    TCR = snakemake.input.tcr
    BRC = snakemake.input.brc
    output = snakemake.output.cat
except:
    parser = get_argparser()
    args = parser.parse_args()
    
    TCR = args.tcr
    BRC = args.barcodes
    output = args.output

# ## Load
tcr_df = pd.read_csv(TCR, low_memory=False)
brc_df = pd.read_csv(BRC, low_memory=False)

# ## Merge clonotypes and barcodes
clonotype_barcode_specificity_df = pd.merge(tcr_df,
                                            brc_df,
                                            how='right', on='gem')

columns = clonotype_barcode_specificity_df.columns.to_list()
columns.remove('umi_count_lst_TRA')
columns.remove('umi_count_lst_TRB')
columns.remove('umi_count_lst_mhc')
columns.remove('umi_count_lst_cd8')
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


# ## Write to excel
clonotype_barcode_specificity_df.to_csv(output, index=False)


