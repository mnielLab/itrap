#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import json
import itertools
import matplotlib.pyplot as plt

from D_plot_specificity_matrix_utils import epitope_sorter_index

plt.style.use('ggplot')


# # Args

#EXP = "exp5"
#PLATFORM = "ILLUMINA"
#MAPPING = 'KMA' # BLAST
BARCODE_SYSTEM = snakemake.params[0] #'10x' #'AKB' #

if BARCODE_SYSTEM == '10x':
    BARCODE_SYSTEM_REGEX = "^(?!.*A\d+B\d+).*$"
    ANTIBODY_REGEX = "HASH"
if BARCODE_SYSTEM == 'AKB':
    BARCODE_SYSTEM_REGEX = "^A\d+B\d+"
    ANTIBODY_REGEX = "A4000"


# ## Input data
map_file = snakemake.input[0]

# OBS! We wont always have response data?!
specificity_annotations = snakemake.input[1] #"/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/barcode_library/barcode_specificity_annotations.tab"
response_annotations = snakemake.input[2] #"/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp3_MHC/barcode_library/detected_responses_annotation.xlsx"
umi_annotations = snakemake.input[3] #"/Volumes/tuba/kamilla/10x-barcoding/results/kma_parser_parallel.tsv"

# ## Output data
output = snakemake.output[0]
report = snakemake.output[1]

# # Import input
map_df = pd.read_csv(map_file) #, usecols=['query_id', 'template_id', 'gem', 'bit_score', 'alignment_length', 'tso', 'b_primer', 'anneal', 'a_primer', 'match'], sep=" ", names=["read_id", "gem", "tso", "b_primer", "anneal", "cd8_primer", "mhc_primer"]
specificity_df = pd.read_excel(specificity_annotations, sheet_name=0, names=['barcode', 'peptide', 'HLA', 'epitope', 'epitope_rank']) #skiprows=1
cd8_specificity_df = pd.read_excel(specificity_annotations, sheet_name=1, names=['barcode','sample_id','HLA_A','HLA_B','HLA_C','comment'])
response_df = pd.read_excel(response_annotations, usecols=['barcode_cd8', 'peptide'])
umi_df = pd.read_csv(umi_annotations, sep='\t', usecols=['read', 'gem', 'A_N6', 'B_N6'])

# Report raw counts
rprt = {'total_barcode_mapped_gems': list(set(map_df.gem)), 'total_umi_mapped_gems': list(set(umi_df.gem + "-1"))}

# Preprocess UMI table
umi_df.rename(columns={'read':'query_id'}, inplace=True)
umi_df.fillna('', inplace=True)
umi_df = umi_df[(umi_df.A_N6.apply(lambda x: len(str(x))==6)) & (umi_df.B_N6.apply(lambda x: len(str(x))==6))]
umi_df['umi'] = umi_df.A_N6 + umi_df.B_N6

rprt['full_length_umi_gems'] = list(set(umi_df.gem + "-1"))

# ## Process mapping table
print("All reads and GEMs")
print("Reads: %i" %map_df.shape[0])
print("GEMs: %i" %map_df.gem.unique().shape[0])
map_df = map_df[(map_df.credible_alignment == True) & (map_df.barcode.str.contains(BARCODE_SYSTEM_REGEX))]
rprt['credibly_mapped_gems'] = list(set(map_df.gem))
print("Credible reads and GEMs")
print("Reads: %i" %map_df.shape[0])
print("GEMs: %i" %map_df.gem.unique().shape[0])
map_df = pd.merge(map_df, umi_df[['query_id', 'umi']], on='query_id', how='inner')
rprt['credibly_mapped_with_umi_gems'] = list(set(map_df.gem))
print("Credible, UMI annotated reads and GEMs")
print("Reads: %i" %map_df.shape[0])
print("GEMs: %i" %map_df.gem.unique().shape[0])

assert map_df.barcode.isna().sum() == 0
assert map_df.umi.isna().sum() == 0

# Preprocess specificity table
specificity_df['peptide'] = specificity_df.peptide.str.strip().str.split("_", expand=True)[0]
specificity_df['peptide_HLA'] = specificity_df.peptide + ' ' + specificity_df.HLA #.str.split('_').str[0]

cd8_specificity_df.fillna('', inplace=True)
cd8_specificity_df['HLA'] = cd8_specificity_df['HLA_A'] + ', ' + cd8_specificity_df['HLA_B'] + ', ' + cd8_specificity_df['HLA_C']
cd8_specificity_df['HLA'] = cd8_specificity_df.HLA.str.split(', ').apply(lambda x: [i for i in x if i!= ''])
CDX_NOMENCLATURE = cd8_specificity_df.barcode
# At this point I have a table with multiple lines per GEM: each line corresponds to the best annotated read. The annotation of reads may agree on the same barcode or may disagree. Later I will count the number of reads for each barcode and only present the barcode with most reads.
# Reads are filtered so that only reads with full length UMI (12 bp) are represented. Thus approximately 18.000 reads were removed of which 17.000 were due to lack of complete UMI.


# ## Functions
def annotate_lst_per_template(var): # umi, query_id
    dct = map_df.groupby(['gem', 'template_id'])[var].unique().to_dict()
    return map_df.set_index(['gem', 'template_id']).index.map(dct)

def annotate_count(df, var): #umi_lst, read_lst
    return df[var].apply(len)

#def annotate_count(var): #umi_lst, read_lst
#    return map_df[var].apply(lambda x: len(x))

def annotate_lst(df, var): # template_id, epitope, peptide, peptide_HLA, HLA
    if all(df.applymap(type)[var] == list):
        dct = df.groupby('gem')[var].apply(list).apply(lambda x: list(k for k,_ in itertools.groupby(x))).to_dict()
    else:
        dct = df.groupby('gem')[var].unique().to_dict()
    return df.gem.map(dct)

#def annotate_template_lst(df):
#    dct = df.groupby(['gem']).template_id.unique().to_dict()
#    return df.gem.map(dct)

def annotate_count_lst(df, var): #umi_count, read_count
    dct = df.drop_duplicates(subset=['gem','template_id']).groupby(['gem'])[var].apply(list).to_dict()
    return df.gem.map(dct)

def annotate_single_barcode(df):
    return df.umi_count_lst.apply(lambda x: True if len(x)==1 else False)

def annotate_peptide_specificities():
    return pd.merge(mhc_df, specificity_df, how='left', on='barcode')

def annotate_sample_specificities():
    return pd.merge(cd8_df, cd8_specificity_df, how='left', on='barcode')

#def annotate_epitope_lst(df):
#    dct = df.groupby(['gem']).epitope.unique().to_dict()
#    return df.gem.map(dct)

#def annotate_specificity_lst(df):
#    dct = df.groupby(['gem']).peptide_HLA.unique().to_dict()
#    return df.gem.map(dct)

#def annotate_multiplets(df):
#	return df.umi_count_lst.apply(len)

def annotate_delta_umi(df):
    #def calc_delta(x):
    #    if len(x) == 1:
    #        return 100
    #    else:
    #        return int((x[-1]-x[-2])/float(x[-1])*100)
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


# ## Annotate barcode UMI and read counts
map_df['umi_lst'] = annotate_lst_per_template('umi')
map_df['umi_count'] = annotate_count(map_df, 'umi_lst') #annotate_count('umi_lst')

map_df['read_lst'] = annotate_lst_per_template('query_id')
map_df['read_count'] = annotate_count(map_df, 'read_lst')

#cd8_df = map_df[map_df.template_id.str.contains(ANTIBODY_REGEX, na = False)].copy()
#mhc_df = map_df[~map_df.template_id.str.contains(ANTIBODY_REGEX, na = False)].copy()

print(map_df.template_id.str.split('_', expand=True)[0].isin(CDX_NOMENCLATURE).any())
cd8_df = map_df[map_df.template_id.str.split('_', expand=True)[0].isin(CDX_NOMENCLATURE)].copy()
mhc_df = map_df[~map_df.template_id.str.split('_', expand=True)[0].isin(CDX_NOMENCLATURE)].copy()

mhc_df.sort_values(by=['gem','umi_count','score','alignment_length'], inplace=True)
cd8_df.sort_values(by=['gem','umi_count','score','alignment_length'], inplace=True)

cd8_df['template_lst'] = annotate_lst(cd8_df, 'template_id') #annotate_template_lst(cd8_df)
mhc_df['template_lst'] = annotate_lst(mhc_df, 'template_id') #annotate_template_lst(mhc_df)

cd8_df['umi_count_lst'] = annotate_count_lst(cd8_df, 'umi_count')
mhc_df['umi_count_lst'] = annotate_count_lst(mhc_df, 'umi_count')

cd8_df['read_count_lst'] = annotate_count_lst(cd8_df, 'read_count')
mhc_df['read_count_lst'] = annotate_count_lst(mhc_df, 'read_count')

cd8_df['single_barcode'] = annotate_single_barcode(cd8_df)
mhc_df['single_barcode'] = annotate_single_barcode(mhc_df)

cd8_df = annotate_sample_specificities()
mhc_df = annotate_peptide_specificities()

cd8_df['HLA_lst'] = annotate_lst(cd8_df, 'HLA')
cd8_df['sample_id_lst'] = annotate_lst(cd8_df, 'sample_id')
cd8_df['HLA_pool_cd8'] = cd8_df.gem.map(cd8_df.groupby('gem')['HLA'].apply(lambda x: np.unique([z for y in x for z in y])).to_dict())

mhc_df['epitope_lst'] = annotate_lst(mhc_df, 'epitope') #annotate_epitope_lst(mhc_df)
mhc_df['peptide_lst'] = annotate_lst(mhc_df, 'peptide') #annotate_peptide_lst(mhc_df)
mhc_df['HLA_lst'] = annotate_lst(mhc_df, 'HLA')
mhc_df['peptide_HLA_lst'] = annotate_lst(mhc_df, 'peptide_HLA') #annotate_specificity_lst(mhc_df)

cd8_df['multiplets'] = annotate_count(cd8_df, 'umi_count_lst') #annotate_multiplets(cd8_df) 
mhc_df['multiplets'] = annotate_count(mhc_df, 'umi_count_lst') #annotate_multiplets(mhc_df)

cd8_df['delta_umi'] = annotate_delta_umi(cd8_df)
mhc_df['delta_umi'] = annotate_delta_umi(mhc_df)

# QC
print("Credible, UMI annotated reads and GEMs")
print("Reads: %i" %map_df.shape[0])
print("GEMs: %i" %map_df.gem.unique().shape[0])

print("MHC barcodes")
print("Reads: %i" %mhc_df.shape[0])
print("GEMs: %i" %mhc_df.gem.unique().shape[0])

print("Antibody barcodes")
print("Reads: %i" %cd8_df.shape[0])
print("GEMs: %i" %cd8_df.gem.unique().shape[0])

assert map_df.shape[0] == (mhc_df.shape[0] + cd8_df.shape[0])

# ## Collapse into 1 Gem per row
unique_cd8_df = cd8_df.drop_duplicates(subset=['gem'], keep='last')
unique_mhc_df = mhc_df.drop_duplicates(subset=['gem'], keep='last')


# QC
print("MHC: 1 GEM / row")
print("Rows: %i" %unique_mhc_df.shape[0])
print("GEMs: %i" %unique_mhc_df.gem.unique().shape[0])

print("Antibody: 1 GEM / row")
print("Rows: %i" %unique_cd8_df.shape[0])
print("GEMs: %i" %unique_cd8_df.gem.unique().shape[0])

rprt['credibly_mapped_with_umi_mhc_gems'] = list(set(unique_mhc_df.gem))
rprt['credibly_mapped_with_umi_cdx_gems'] = list(set(unique_cd8_df.gem))


# ## MERGE
barcode_df = pd.merge(unique_mhc_df, unique_cd8_df, how='outer', on='gem', suffixes=('_mhc', '_cd8'))

barcode_df['detected_response'] = annotate_detected_response(barcode_df)
barcode_df['peptide_assayed'] = annotate_peptide_assayed(barcode_df)

# QC
print("Merge of MHC and antibody barcodes")
print("Rows: %i" %barcode_df.shape[0])
print("GEMs: %i" %barcode_df.gem.unique().shape[0])
assert map_df.gem.unique().shape[0] == barcode_df.gem.unique().shape[0]

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

# ## Check that annotated peptide HLA matches CDX HLA annotation
barcode_df['HLA_match'] = barcode_df.apply(lambda row: row.HLA_mhc in row.HLA_cd8 if (row.HLA_mhc==row.HLA_mhc) & (type(row.HLA_cd8) == list) else np.nan, axis=1)
barcode_df['likely_HLA_mhc'] = barcode_df.apply(lambda row: get_likely_targets(row), axis=1)

# ## Specificity matrix
specificity_matrix = unique_mhc_df.dropna(subset=['epitope']).pivot(index='gem', columns='epitope', values='umi_count')

# ## MERGE
barcode_specificity_df = pd.merge(barcode_df, #response_df, how='left', on=['barcode_cd8', 'peptide']).merge(
								  specificity_matrix, how='left', on='gem')


#barcode_specificity_df['detected_response'] = annotate_detected_response(barcode_specificity_df)
#barcode_specificity_df['peptide_assayed'] = annotate_peptide_assayed(barcode_specificity_df)

# Check epitope ranking
if not any(barcode_specificity_df.epitope_rank):
    barcode_specificity_df['epitope_rank'] = epitope_sorter_index(barcode_specificity_df)

print("Merge of barcodes and specificity matrix")
print("Rows: %i" %barcode_specificity_df.shape[0])
print("GEMs: %i" %barcode_specificity_df.gem.unique().shape[0])
assert map_df.gem.unique().shape[0] == barcode_specificity_df.gem.unique().shape[0]

# ## Write data
new_column_order = ['gem',
 'template_id_mhc',
 'template_lst_mhc',
 'single_barcode_mhc',
 'umi_count_mhc',
 'umi_count_lst_mhc',
 'read_count_mhc',
 'read_count_lst_mhc',
 'multiplets_mhc',
 'delta_umi_mhc',
 'template_id_cd8',
 'template_lst_cd8',
 'single_barcode_cd8',
 'umi_count_cd8',
 'umi_count_lst_cd8',
 'read_count_cd8',
 'read_count_lst_cd8',
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
 'epitope_rank'] #+ specificity_matrix.columns.to_list()

barcode_specificity_df[new_column_order].to_csv(output, index=False)

with open(report, 'w') as outfile:
    json.dump(rprt, outfile)





