#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
#import itertools
from ast import literal_eval
import re
#from random import sample
#import random
import re
import yaml

#plt.style.use('ggplot')

def HLA_cd8_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace(",", "").replace("'","").split(" ")

def cdr3_lst_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace("'","").split(" ")

def epitope_converter(x):
    #define format of datetime
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
    #tmp = df[idx1 & idx2]
    dct = df.groupby(['ct','peptide_HLA']).gem.count() > 1
    idx = df.set_index(['ct','peptide_HLA']).index.map(dct)
    return idx.fillna(False)

##########################################################
#                         Inputs                         #
##########################################################
VALID = snakemake.input.valid
THRESHOLD = snakemake.input.opt_thr
GEX = snakemake.input.gex #'../experiments/exp13/run1_archive/tcr/usable_gems.txt'

filter_set = snakemake.params.flt
##########################################################
#                         Output                         #
##########################################################
OUTPUT = snakemake.output.flt

##########################################################
#                          Load                          #
##########################################################

opt_thr = pd.read_csv(THRESHOLD, index_col=0, header=None, names=['thr']).thr.dropna()
gex = pd.read_csv(GEX, header=None, names=['gem'])
df = pd.read_csv(VALID, converters=converters, low_memory=False)


##########################################################
#                          Prep                          #
##########################################################
# when filtering for optimal threshold it is important to have values in UMI and delta
df.fillna({'umi_count_mhc':0, 'delta_umi_mhc':0, "umi_count_mhc_rel":0,
           'umi_count_TRA':0, 'delta_umi_TRA':0,
           'umi_count_TRB':0, 'delta_umi_TRB':0}, inplace=True)

# Add extra features
df['gex'] = df.gem.isin(gex.gem)
df.single_barcode_mhc = np.where(df.single_barcode_mhc, 'pMHC singlet','pMHC multiplet')
df['clonotype_multiplet'] = df.ct.map(df.groupby('ct').size() > 1)
df['HLA_match_per_gem'] = df.apply(lambda row: row.HLA_mhc in row.HLA_cd8 if row.HLA_cd8 == row.HLA_cd8 else False, axis=1)


##########################################################
#                         Filters                        #
##########################################################
idx0 = ~df.gem.isna() # Total
idx1 = eval(' & '.join([f'(df.{k} >= {abs(v)})' for k,v in opt_thr.items()])) # optimal threshold
idx2 = df.hto_global_class == 'Singlet' # Hashing singlets
idx3 = df.apply(lambda row: row.peptide_HLA.split()[-1] in row.HLA_cd8 if (notnan(row.peptide_HLA) & notnan(row.HLA_cd8)) else False, axis=1) # Matching HLA
idx4 = df['exclude_single-chain_TCRs'] # Complete TCRs
idx5 = get_multiplets(df) # Only specificities in multiplets
try:
    idx6 = df.cell_flag # is_cell
except AttributeError:
    idx6 = df.gem.isna() # All false
#idx7 = df.is_cell_gex
idx8 = df.gex # Gene expression data

is_cell_gex = ' (GEX)' if any(idx8) else ''

if filter_set == 'indv':
    # Showing individual effects of filtering
    filterings = [idx0,
              idx1,
              idx3,
              idx2,
              idx4,
              idx5,
              idx6,
              #idx7,
              idx8]
    labels = ['total','optimal threshold',
          'matching HLA',
          'hashing singlets',
          'complete TCRs',
          'specificity multiplets',
          'is cell%s' %is_cell_gex,
          'is viable cell']
    palette = ['grey','yellow','#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','black']
    
    flt_to_remove = list()
    for i, flt in enumerate(filterings):
        print(labels[i], sum(flt))
        if sum(flt) == 0:
            print('removed:', labels[i])
            flt_to_remove.append(i)
            
    for i in flt_to_remove[::-1]:
        del filterings[i]
        del labels[i]
        del palette[i]

elif filter_set == 'comb':
    # Showing combined effects in the same order
    #filter_set = 'comb'
    labels = ['total','optimal threshold',
              'matching HLA',
              'hashing singlets',
              'complete TCRs',
              'specificity multiplets',
              'is cell%s' %is_cell_gex,
              'is viable cell']
    palette = ['grey','yellow','#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','black'] #'#253494'
    
    flt_to_remove = list()
    filterings = [idx0]
    for i, flt in enumerate([idx1, idx3, idx2, idx4, idx5, idx6, idx8], start=1):
        remaining_gems = sum(filterings[-1] & flt)
        print(labels[i], remaining_gems)
        if remaining_gems > 0:
            filterings.append((filterings[-1] & flt))
        else:
            print('removed:', labels[i])
            flt_to_remove.append(i)
            
    for i in flt_to_remove[::-1]:
        del labels[i]
        del palette[i]
            
    #filterings = [idx0,
    #              idx1,
    #              (idx1 & idx3),
    #              (idx1 & idx2 & idx3),
    #              (idx1 & idx2 & idx3 & idx4),
    #              (idx1 & idx2 & idx3 & idx4 & idx5),
    #              (idx1 & idx2 & idx3 & idx4 & idx5 & idx6),
    #              #(idx1 & idx2 & idx3 & idx4 & idx5 & idx7),
    #              (idx1 & idx2 & idx3 & idx4 & idx5 & idx8)]
    
else:
    print('filterset name unknown')

##########################################################
#                       Output prep                      #
##########################################################
dct = dict(labels = labels,
           platte = palette)

tmp = pd.concat(filterings)

##########################################################
#                  Write output to file                  #
##########################################################
with open(OUTPUT, 'w') as outfile:
    yaml.dump(dct, outfile)
