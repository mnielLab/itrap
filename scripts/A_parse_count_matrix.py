#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from scipy.io import mmread
import itertools
import argparse

def get_argparser():
    """Return an argparse argument parser."""
    parser = argparse.ArgumentParser(prog = 'Clean & reformat barcode data',
                                     description = 'Reformats Cellranger barcode output')
    add_arguments(parser)
    return parser

def add_arguments(parser):
    parser.add_argument('--features', required=True, nargs=1, help='Cellranger multi output: features.tsv.gz')
    parser.add_argument('--barcodes', required=True, nargs=1, help='Cellranger multi output: barcodes.tsv.gz')
    parser.add_argument('--matrix', required=True, nargs=1, help='Cellranger multi output: matrix.tsv.gz')
    parser.add_argument('--multimers', required=True, help='Filepath for multimer barcode specifications')
    parser.add_argument('--responses', required=False, default=None, help='Filepath for previously detected responses')
    parser.add_argument('--labels', required=False, nargs=1, default=None, help='Filename for GEMs labelled by FACS sorting')
    parser.add_argument('--output', required=True, help='Filepath of output data')
    parser.add_argument('--mhc_custom', action='store_true', help='Flag if pMHC barcodes are of custom design (non-immundex)')
    parser.add_argument('--hsh_custom', action='store_true', help='Flag if cell hashing barcodes are of custom design (non-immundex)')
    parser.add_argument('--mrk_custom', action='store_true', help='Flag if surface marker barcodes are of custom design (non-immundex)')

def annotate_count(df, var):
    return df[var].apply(len)

def annotate_lst(df, var):
    if all(df.applymap(type)[var] == list):
        dct = df.groupby('gem')[var].apply(list)
    else:
        dct = df.groupby('gem')[var].apply(list)
    return df.gem.map(dct)

def annotate_pool(df):
    dct = df.groupby('gem').HLA.apply(lambda x: np.unique([z for y in x for z in y]))
    return df.gem.map(dct)

def annotate_delta_umi(df):
    def calc_delta(x):
        if len(x) == 1:
            return x[-1]/0.25
        elif len(x) == 0:
            return 0
        else:
            return (x[-1])/(x[-2]+0.25)
    return df.umi_count_lst.apply(calc_delta)

def annotate_detected_response(df):
    dct = response_df.groupby(['peptide','barcode_cd8']).apply(any).to_dict()
    return df.set_index(['peptide','barcode_cd8']).index.map(dct)

def annotate_peptide_assayed(df):
    return np.where(df.peptide.isin(response_df.peptide), True, False)

def get_likely_targets(row):
    from itertools import compress
    if (type(row.HLA_cd8) is list) & (type(row.HLA_lst_mhc) is list):
        chec = [item in row.HLA_cd8 for item in row.HLA_lst_mhc]
        idxs = list(compress(range(len(chec)), chec))
        if idxs == []:
            return np.nan
        else:
            return [row.HLA_lst_mhc[i] for i in idxs]
    else:
        return np.nan


try:
    BRC = snakemake.input.brc 
    GEM = snakemake.input.gem
    MTX = snakemake.input.mtx 
    BARCODE_SOURCE = snakemake.input.cnt
    BARCODE_ANNOTATIONS = snakemake.input.lbl
    RESPONSE_ANNOTATIONS = snakemake.input.ann
    MHC_CUSTOM = snakemake.params.mhc_custom
    HSH_CUSTOM = snakemake.params.hsh_custom
    MRK_CUSTOM = snakemake.params.mrk_custom
    OUTPUT = snakemake.output[0]
except:
    parser = get_argparser()
    args = parser.parse_args()
    
    BRC = args.features
    GEM = args.barcodes
    MTX = args.matrix
    BARCODE_SOURCE = args.labels
    BARCODE_ANNOTATIONS = args.multimers
    RESPONSE_ANNOTATIONS = args.responses
    MHC_CUSTOM = args.mhc_custom
    HSH_CUSTOM = args.hsh_custom
    MRK_CUSTOM = args.mrk_custom
    OUTPUT = args.output

# # Load
response_df = pd.read_excel(RESPONSE_ANNOTATIONS, usecols=['barcode_cd8', 'peptide'])

sources = dict()
if BARCODE_SOURCE:
    for filename in BARCODE_SOURCE:
        tmp_df = pd.read_csv(filename, header=None, names=['gem','barcode','umi','label'])
        tmp_df.drop_duplicates(subset=['gem','barcode'], keep='first', inplace=True) # why is this necessay?
        if '10x' in filename:
            tmp_df.set_index(['gem'], inplace=True)
            sources['10x'] = tmp_df.label
        else:
            tmp_df.set_index(['gem','barcode'], inplace=True)
            sources['kma'] = tmp_df.label

labels = dict()
barcodes = list()
for label, sheet in zip(['mhc','hsh','mrk'],['MHC','HSH','MRK']):
    labels[label] = pd.read_excel(BARCODE_ANNOTATIONS, sheet_name=sheet)
    
    if label == 'mhc':
        labels[label]['peptide'] = labels[label].peptide.str.strip().str.split("_", expand=True)[0]
        labels[label]['peptide_HLA'] = labels[label].peptide + ' ' + labels[label].HLA
        barcodes += labels[label]['barcode'].astype(str).to_list()
        
    if label == 'hsh':
        labels[label].fillna('', inplace=True)
        labels[label]['HLA'] = labels[label]['HLA_A'] + ', ' + labels[label]['HLA_B'] + ', ' + labels[label]['HLA_C']
        labels[label]['HLA'] = labels[label]['HLA'].str.split(r',\s?').apply(lambda x: [i for i in x if i!= ''])
        barcodes += labels[label]['barcode'].astype(str).to_list() 
        
    if label == 'mrk':
        continue

t = list()
for B,G,M in zip(BRC, GEM, MTX):
    g = np.loadtxt(G, dtype='U36')
    b = pd.read_csv(B, sep='\t', header=None, names=['barcode','name','feature'])
    m = mmread(M)
    m = pd.DataFrame.sparse.from_spmatrix(m, index=b.barcode, columns=g)
    m = m.loc[m.index.isin(barcodes),:]
    t.append(m.sparse.to_dense().T)

wide_df = pd.concat(t, axis=1) # Concatenating along rows: GEMs

for d in t:
    del d
del t

# # Prep
wide_df.replace(0, np.nan, inplace=True)

long_df = wide_df.melt(value_vars=wide_df.columns, ignore_index=False)
long_df.reset_index(inplace=True)
long_df.rename(columns={'index':'gem','value':'umi_count'}, inplace=True)
long_df.dropna(inplace=True)

dfs = dict()

long_df.barcode = long_df.barcode.astype(str)
for key, frame in labels.items():
    frame.barcode = frame.barcode.astype(str)
    
    df = pd.merge(long_df, frame, on='barcode')
    df.sort_values(by=['gem','umi_count'], inplace=True)
    
    df['umi_count_lst'] = annotate_lst(df, 'umi_count')
    df['delta_umi'] = annotate_delta_umi(df)
    df['brc_count'] = df.umi_count_lst.apply(len)
    df['single_barcode'] = df.brc_count == 1
    df['multiplets'] = df.brc_count > 1
    df['template_lst'] = annotate_lst(df, 'barcode')
    
    if key == 'mhc':
        df['epitope_lst'] = annotate_lst(df, 'epitope')
        df['peptide_lst'] = annotate_lst(df, 'peptide')
        df['peptide_HLA_lst'] = annotate_lst(df, 'peptide_HLA')
        df['HLA_lst'] = annotate_lst(df, 'HLA')
        df['HLA_pool'] = annotate_pool(df) # Not really necessary for mhc
        
        if MHC_CUSTOM:
            df['label'] = df.set_index(['gem','barcode']).index.map(sources['kma'])
        else:
            if sources:
                df['label'] = df.gem.map(sources['10x'])
            else:
                df['label'] = 1
        
    if key == 'hsh':
        df['sample_id_lst'] = annotate_lst(df, 'sample_id')
        df['HLA_lst'] = annotate_lst(df, 'HLA')
        df['HLA_pool'] = annotate_pool(df)
        
        assert all(df.sample_id_lst.apply(len) == df.HLA_lst.apply(len))
        assert all(df.HLA_lst.apply(len) == df.umi_count_lst.apply(len))
        
        if HSH_CUSTOM:
            df['label'] = df.set_index(['gem','barcode']).index.map(sources['kma'])
        else:
            if sources:
                df['label'] = df.gem.map(sources['10x'])
            else:
                df['label'] = 1
        
    if key == 'mrk':
        df['marker_lst'] = annotate_lst(df, 'marker')
        
        if MRK_CUSTOM:
            df['label'] = df.set_index(['gem','barcode']).index.map(sources['kma'])
        else:
            if sources:
                df['label'] = df.gem.map(sources['10x'])
            else:
                df['label'] = 1
            
        df.columns = [name if name in ['gem','marker','marker_lst'] else name + '_mrk' for name in df.columns] # Hack
        
    df.drop_duplicates(subset=['gem'], keep='last', inplace=True)
    if key == 'hsh':
        assert all(df.sample_id_lst.apply(len) == df.HLA_lst.apply(len))
        assert all(df.HLA_lst.apply(len) == df.umi_count_lst.apply(len))
    dfs[key] = df
    
df = (dfs['mhc']
      .merge(dfs['hsh'], on='gem', how='outer', suffixes=('_mhc','_cd8'))
      .merge(dfs['mrk'], on='gem', how='outer', suffixes=('','_mrk')))
    
# ## Check that annotated peptide HLA matches CDX HLA annotation
df['HLA_match'] = df.apply(lambda row: row.HLA_mhc in row.HLA_cd8 if (row.HLA_mhc==row.HLA_mhc) & (type(row.HLA_cd8) == list) else np.nan, axis=1)
df['likely_HLA_mhc'] = df.apply(lambda row: get_likely_targets(row), axis=1)

df['detected_response'] = annotate_detected_response(df)
df['peptide_assayed'] = annotate_peptide_assayed(df)

# ## Write data
new_column_order = ['gem',
 'umi_count_mhc',
 'umi_count_lst_mhc',
 'delta_umi_mhc',
 'single_barcode_mhc',
 'multiplets_mhc',
 'umi_count_cd8',
 'umi_count_lst_cd8',
 'delta_umi_cd8',
 'single_barcode_cd8',
 'multiplets_cd8',
 'detected_response',
 'peptide_assayed',
 'sample_id',
 'sample_id_lst',
 'HLA_pool_cd8',
 'HLA_lst_cd8',
 'HLA_cd8',
 'HLA_match',
 'HLA_mhc',
 'HLA_lst_mhc',
 'likely_HLA_mhc',
 'peptide',
 'peptide_lst',
 'peptide_HLA',
 'peptide_HLA_lst',
 'epitope',
 'epitope_lst',
 'epitope_rank',
 'marker',
 'marker_lst',
 'umi_count_lst_mrk',
 'delta_umi_mrk',
 'single_barcode_mrk',
 'multiplets_mrk',
 'label_mhc','label_cd8','label_mrk']

df[new_column_order].to_csv(OUTPUT, index=False)







