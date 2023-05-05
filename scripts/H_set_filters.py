#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from ast import literal_eval
import re
import yaml
import argparse

def get_argparser():
    """Return an argparse argument parser."""
    parser = argparse.ArgumentParser(prog = 'Set filters',
                                     description = 'Generates boolean arrays by which to filter data on.')
    add_arguments(parser)
    return parser

def add_arguments(parser):
    parser.add_argument('--data', required=True, help='Filepath for data')
    parser.add_argument('--opt-thr', required=True, help='Filepaths for optimal threshold')
    parser.add_argument('--setting', required=True, choices=['indv','comb'], help='Defined whether filters are applied individually or combined additively.')
    parser.add_argument('--labels', required=True, help='Filepath for output filter labels')
    parser.add_argument('--filters', required=True, help='Filepath for output filters')

def HLA_cd8_converter(x):
    return x.replace("[","").replace("]","").replace(",", "").replace("'","").split(" ")

def cdr3_lst_converter(x):
    return x.replace("[","").replace("]","").replace("'","").split(" ")

def epitope_converter(x):
    return [y for y in x.replace("[","").replace("]","").replace("\n","").split("'") if (y != '') & (y != ' ')]

def peptide_hla_converter(x):
    return re.findall("\w+\s{1}\w{1}\d+", x.replace("[","").replace("]","").replace("\n","").replace("'",""))

def literal_converter(val):
    # replace NaN with '' and perform literal eval on the rest
    return [] if val == '' else literal_eval(val)

converters = {'peptide_HLA_lst': peptide_hla_converter,
              'umi_count_lst_mhc': literal_eval,
              'umi_count_lst_TRA': literal_converter,'umi_count_lst_TRB': literal_converter,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'cdr3_lst_TRB': cdr3_lst_converter,
              'HLA_lst_mhc': cdr3_lst_converter,'HLA_cd8': HLA_cd8_converter} #

def notnan(x):
    return x == x

def get_multiplets(df):
    dct = df.groupby(['ct','peptide_HLA']).gem.count() > 1
    idx = df.set_index(['ct','peptide_HLA']).index.map(dct)
    return pd.Series(idx.fillna(False))

##########################################################
#                          Load                          #
##########################################################

try:
    VALID = snakemake.input.valid
    THRESHOLD = snakemake.input.opt_thr
    filter_set = snakemake.params.flt
    YAML = snakemake.output.lbl
    DATA = snakemake.output.flt
except:
    parser = get_argparser()
    args = parser.parse_args()
    
    VALID = args.data
    THRESHOLD = args.opt_thr
    filter_set = args.setting
    YAML = args.labels
    DATA = args.filters
    

opt_thr = pd.read_csv(THRESHOLD, index_col=0, header=None, names=['thr']).thr.dropna()
df = pd.read_csv(VALID, converters=converters, low_memory=False)


##########################################################
#                          Prep                          #
##########################################################
# when filtering for optimal threshold it is important to have values in UMI and delta
df.fillna({'umi_count_mhc':0, 'delta_umi_mhc':0, "umi_count_mhc_rel":0,
           'umi_count_TRA':0, 'delta_umi_TRA':0,
           'umi_count_TRB':0, 'delta_umi_TRB':0}, inplace=True)

# Add extra features
df.single_barcode_mhc = np.where(df.peptide_HLA_lst.apply(len) > 1, 'pMHC singlet','pMHC multiplet')
df['clonotype_multiplet'] = df.ct.map(df.groupby('ct').size() > 1)
df['HLA_match_per_gem'] = df.apply(lambda row: row.HLA_mhc in row.HLA_cd8 if row.HLA_cd8 == row.HLA_cd8 else False, axis=1)
df['complete_tcrs'] = df.cdr3_TRA.notna() & df.cdr3_TRB.notna()


##########################################################
#                         Filters                        #
##########################################################
# idx0: raw
# idx1: UMI thresholds
# idx2: Hashing singlets
# idx3: Matching HLA
# idx4: Complete TCRs
# idx5: Specificity multiplets
# idx6: Is cell (Cellranger)
# idx7: Viable cells (GEX)
idx0 = ~df.gem.isna()
idx1 = eval(' & '.join([f'(df.{k} >= {abs(v)})' for k,v in opt_thr.items()]))
idx2 = df.hto_global_class == 'Singlet'
idx3 = df.apply(lambda row: row.peptide_HLA.split()[-1] in row.HLA_cd8 if (notnan(row.peptide_HLA) & notnan(row.HLA_cd8)) else False, axis=1)
idx4 = df['complete_tcrs'] #exclude_single-chain_TCRs
idx5 = get_multiplets(df)
try:
    idx6 = df.cell_flag # is_cell
except AttributeError:
    idx6 = df.gem.isna() # All false
idx7 = df.gex

if filter_set == 'indv':
    # Showing individual effects of filtering
    filterings = [idx0,
              idx1,
              idx3,
              idx2,
              idx4,
              idx5,
              idx6,
              idx7]
    labels = ['total','optimal threshold',
          'matching HLA',
          'hashing singlets',
          'complete TCRs',
          'specificity multiplets',
          'is cell%s' %(' (GEX)' if any(idx6) else ''),
          'is viable cell']
    palette = ['grey','yellow','#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','black']
    
    flt_to_remove = list()
    for i, flt in enumerate(filterings):
        if sum(flt) == 0:
            flt_to_remove.append(i)
            
    for i in flt_to_remove[::-1]:
        del filterings[i]
        del labels[i]
        del palette[i]

elif filter_set == 'comb':
    # Showing combined effects in the same order
    labels = ['total','optimal threshold',
              'matching HLA',
              'hashing singlets',
              'complete TCRs',
              'specificity multiplets',
              'is cell%s' %(' (GEX)' if any(idx6) else ''),
              'is viable cell']
    palette = ['grey','yellow','#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','black'] #'#253494'
    
    flt_to_remove = list()
    filterings = [idx0]
    for i, flt in enumerate([idx1, idx3, idx2, idx4, idx5, idx6, idx7], start=1):
        remaining_gems = sum(filterings[-1] & flt)
        if remaining_gems > 0:
            filterings.append((filterings[-1] & flt))
        else:
            flt_to_remove.append(i)
            
    for i in flt_to_remove[::-1]:
        del labels[i]
        del palette[i]
    
else:
    print('filterset name unknown')

##########################################################
#                       Output prep                      #
##########################################################
dct = dict(labels = labels,
           palette = palette)

tmp = pd.concat(filterings, axis=1)
tmp.columns = labels

##########################################################
#                  Write output to file                  #
##########################################################
with open(YAML, 'w') as outfile:
    yaml.dump(dct, outfile)
    
tmp.to_csv(DATA, index=False)
