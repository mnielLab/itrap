#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from scipy.io import mmread
import itertools
#from functools import reduce

# ## Functions
def annotate_count(df, var): #umi_lst, read_lst
    return df[var].apply(len)

def annotate_lst(df, var): # template_id, epitope, peptide, peptide_HLA, HLA
    if all(df.applymap(type)[var] == list):
        #dct = df.groupby('gem')[var].apply(list).apply(lambda k: list(k for k,_ in itertools.groupby(k))).to_dict()
        dct = df.groupby('gem')[var].apply(list) #.apply(lambda x: list(k for k,_ in itertools.groupby(x)))
    else:
        dct = df.groupby('gem')[var].apply(list) #.unique()
    return df.gem.map(dct)

def annotate_pool(df):
    # var = HLA
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


# # Input
BRC = snakemake.input.brc 
#'/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run1/mhc/count/features.tsv'
#'/Users/herpov/Documents/mnt/tuba_home/tcr-pmhc-sc-project/experiments/exp13/run2/tcr/cellranger_tot/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz'
GEM = snakemake.input.gem
#'/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run1/mhc/count/barcodes.tsv'
#'/Users/herpov/Documents/mnt/tuba_home/tcr-pmhc-sc-project/experiments/exp13/run2/tcr/cellranger_tot/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz'
MTX = snakemake.input.mtx 
#'/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run1/mhc/count/matrix.mtx'
#'/Users/herpov/Documents/mnt/tuba_home/tcr-pmhc-sc-project/experiments/exp13/run2/tcr/cellranger_tot/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz'
BARCODE_SOURCE = snakemake.input.cnt
BARCODE_ANNOTATIONS = snakemake.input.lbl
#'/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run1/lib/barcode_specificity_annotations.xlsx'
#'/Users/herpov/Documents/mnt/tuba_home/tcr-pmhc-sc-project/experiments/exp13/run2/lib/barcode_specificity_annotations.xlsx'
RESPONSE_ANNOTATIONS = snakemake.input.ann 
#'/home/tuba/herpov/tcr-pmhc-sc-project/tools/detected_responses_annotation.xlsx' 
#'/Users/herpov/Documents/mnt/tuba_home/tcr-pmhc-sc-project/experiments/exp13/run2/lib/detected_responses_annotation.xlsx'

# # Output
OUTPUT = snakemake.output[0] #'tester.csv'

# # Load
response_df = pd.read_excel(RESPONSE_ANNOTATIONS, usecols=['barcode_cd8', 'peptide'])

sources = dict()
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
        print(labels[label])
        labels[label]['HLA'] = labels[label]['HLA_A'] + ', ' + labels[label]['HLA_B'] + ', ' + labels[label]['HLA_C']
        labels[label]['HLA'] = labels[label]['HLA'].str.split(r',\s?').apply(lambda x: [i for i in x if i!= ''])
        barcodes += labels[label]['barcode'].astype(str).to_list() 
        
    if label == 'mrk':
        continue

print(barcodes)

t = list()
for B,G,M in zip(BRC, GEM, MTX):
    g = np.loadtxt(G, dtype='U36')
    b = pd.read_csv(B, sep='\t', header=None, names=['barcode','name','feature'])
    m = mmread(M)
    #m = m.toarray()[b.barcode.isin(barcodes),:]
    m = pd.DataFrame.sparse.from_spmatrix(m, index=b.barcode, columns=g)
    print(m.shape)
    m = m.loc[m.index.isin(barcodes),:]
    print(m.shape)
    t.append(m.sparse.to_dense().T)
    #t.append(pd.DataFrame(m, index=barcodes, columns=g).T) #loc[d.index.isin(barcodes),:].T) #.copy()

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
#long_df['label'] = long_df.gem.map(source_df.label)

dfs = dict()

long_df.barcode = long_df.barcode.astype(str)
for key, frame in labels.items():
    print(key)
    frame.barcode = frame.barcode.astype(str)
    
    df = pd.merge(long_df, frame, on='barcode')
    df.sort_values(by=['gem','umi_count'], inplace=True)
    
    df['umi_count_lst'] = annotate_lst(df, 'umi_count') #umi_lst    
    df['delta_umi'] = annotate_delta_umi(df)
    df['brc_count'] = df.umi_count_lst.apply(len) #brc_count  
    df['single_barcode'] = df.brc_count == 1 #brc_singlet
    df['multiplets'] = df.brc_count > 1 #brc_multiplet
    df['template_lst'] = annotate_lst(df, 'barcode') # barcode == template_id
    
    if key == 'mhc':
        df['epitope_lst'] = annotate_lst(df, 'epitope')
        df['peptide_lst'] = annotate_lst(df, 'peptide')
        df['peptide_HLA_lst'] = annotate_lst(df, 'peptide_HLA')
        df['HLA_lst'] = annotate_lst(df, 'HLA')
        df['HLA_pool'] = annotate_pool(df) # Not really necessary for mhc
        
        if snakemake.params.mhc_custom:
            df['label'] = df.set_index(['gem','barcode']).index.map(sources['kma'])
        else:
            print(df.gem.shape, len(df.gem.unique()))
            print(sources['10x'].shape, len(sources['10x'].index.unique()), sources['10x'].index.duplicated(keep='first').sum())
            print()
            df['label'] = df.gem.map(sources['10x'])
        
    if key == 'hsh':
        df['sample_id_lst'] = annotate_lst(df, 'sample_id')
        df['HLA_lst'] = annotate_lst(df, 'HLA')
        df['HLA_pool'] = annotate_pool(df)
        print(df[['sample_id_lst','HLA_lst','umi_count_lst']])
        assert all(df.sample_id_lst.apply(len) == df.HLA_lst.apply(len))
        assert all(df.HLA_lst.apply(len) == df.umi_count_lst.apply(len))
        
        if snakemake.params.hsh_custom:
            df['label'] = df.set_index(['gem','barcode']).index.map(sources['kma'])
        else:
            df['label'] = df.gem.map(sources['10x'])
        
    if key == 'mrk':
        df['marker_lst'] = annotate_lst(df, 'marker')
        
        if snakemake.params.mrk_custom:
            df['label'] = df.set_index(['gem','barcode']).index.map(sources['kma'])
        else:
            df['label'] = df.gem.map(sources['10x'])
            
        df.columns = [name if name in ['gem','marker','marker_lst'] else name + '_mrk' for name in df.columns] # Hack
        
    df.drop_duplicates(subset=['gem'], keep='last', inplace=True)
    if key == 'hsh':
        assert all(df.sample_id_lst.apply(len) == df.HLA_lst.apply(len))
        assert all(df.HLA_lst.apply(len) == df.umi_count_lst.apply(len))
    dfs[key] = df
    #dfs = dfs.merge(df, on='gem', how='outer', suffixes=('','_%s' % ('cd8' if key == 'hsh' else key)))
    
df = (dfs['mhc']
      .merge(dfs['hsh'], on='gem', how='outer', suffixes=('_mhc','_cd8')) #_hsh
      .merge(dfs['mrk'], on='gem', how='outer', suffixes=('','_mrk')))
    
# ## Check that annotated peptide HLA matches CDX HLA annotation
df['HLA_match'] = df.apply(lambda row: row.HLA_mhc in row.HLA_cd8 if (row.HLA_mhc==row.HLA_mhc) & (type(row.HLA_cd8) == list) else np.nan, axis=1)
df['likely_HLA_mhc'] = df.apply(lambda row: get_likely_targets(row), axis=1)

df['detected_response'] = annotate_detected_response(df)
df['peptide_assayed'] = annotate_peptide_assayed(df)

# ## Write data
new_column_order = ['gem',
 #'template_id_mhc',
 #'template_lst_mhc',
 'umi_count_mhc',
 'umi_count_lst_mhc',
 'delta_umi_mhc',
 #'read_count_mhc',
 #'read_count_lst_mhc',
 'single_barcode_mhc',
 'multiplets_mhc',
 #'template_id_cd8',
 #'template_lst_cd8',
 'umi_count_cd8',
 'umi_count_lst_cd8',
 'delta_umi_cd8',
 'single_barcode_cd8',
 'multiplets_cd8',
 #'read_count_cd8',
 #'read_count_lst_cd8',
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
 'label_mhc','label_cd8','label_mrk'] #+ specificity_matrix.columns.to_list()

df[new_column_order].to_csv(OUTPUT, index=False)







