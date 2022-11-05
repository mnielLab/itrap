#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import yaml
import itertools

import sys
import os
import re

"""
Construct cellranger multi feature file
"""

# # Input files
BARCODE_ANNOTATIONS = snakemake.input.barcode_annotations
BARCODE_LIBRARY = snakemake.input.barcode_library
BRC_INFO = snakemake.input.brc_info # load all_oligos.csv
FQ_PATHS = snakemake.input.fq_paths

path, module = os.path.split(FQ_PATHS)
sys.path.append(path)
import module

# Params imported with module
# pos_tcr_fq_dir
# neg_tcr_fq_dir
# pos_mhc_fq_dir
# neg_mhc_fq_dir
# pos_hsh_fq_dir
# neg_hsh_fq_dir
# pos_mrk_fq_dir
# neg_mrk_fq_dir

# mhc_custom
# hsh_custom
# mrk_custom

# ## Output files
MULTI_CONFIG = snakemake.output.multi_config #"../data/" + EXP + "/barcode_library/barcode_templates.fa"
FEATURE_REFERENCE = snakemake.output.feature_reference

if ('pos' in MULTI_CONFIG) & ('pos' in FEATURE_REFERENCE):
    globals().update(sorting['pos'])
elif ('neg' in MULTI_CONFIG) & ('neg' in FEATURE_REFERENCE):
    globals().update(sorting['neg'])
else:
    print('ERROR')

# Load
mhc_lib = pd.read_csv(BRC_INFO)
hsh_lib = pd.read_excel(BARCODE_LIBRARY, sheet_name=0)
mrk_lib = pd.read_excel(BARCODE_LIBRARY, sheet_name=1)

mhc_ann = pd.read_excel(BARCODE_ANNOTATIONS, sheet_name=0)
hsh_ann = pd.read_excel(BARCODE_ANNOTATIONS, sheet_name=1)
mrk_ann = pd.read_excel(BARCODE_ANNOTATIONS, sheet_name=2)

#bbb = load_barcode_building_blocks()  # No inconvenient

# Write multi config
multi_config = {'gex':f'reference,{snakemake.config.gex_ref}\n',
                'feature':f'reference,{FEATURE_REFERENCE}\n',
                'vdj':f'reference,{snakemake.config.vdj_ref}\n'}

with open(MULTI_CONFIG, 'w') as fh:
    fh.write('[gene-expression]\n')
    if gex:
        fh.write(multi_config['gex'])
    fh.write('[feature]\n')
    if any([mhc_custom==False, hsh_custom==False, mrk_custom==False]):
        fh.write(multi_config['feature'])
    fh.write('[vdj]\n')
    if tcr:
        fh.write(multi_config['vdj'])
    fh.write('[libraries]\n')
    fh.write('fastq_id,fastqs,lanes,feature_types,subsample_rate\n')
    if gex:
        fh.write(f'{gex},{dir_fqs},any,Gene Expression,')
    if tcr:
        fh.write(f'{pos_tcr},{dir_fqs},any,VDJ-T,')
    if mhc & (not mhc_custom):
        fh.write(f'{pos_mhc},{dir_fqs},any,Antigen Capture,')
    if hsh & (not shs_custom):
        fh.write(f'{pos_shs},{dir_fqs},any,Antibody Capture,')
    if mrk & (not mrk_custom):
        fh.write(f'{pos_mrk},{dir_fqs},any,Antibody Capture,')

# Write Feature Reference csv
df_lst = list()
if mhc_custom == False:
    identifier, sequence = 'id_15mer_2','seq_15mer_fw'

    df_tmp = mhc_lib.loc[mhc_lib[identifier].isin(mhc_ann.barcode), [identifier, sequence]]
    df_tmp.columns = ['id','sequence']
    df_tmp['name'] = df_tmp['id'].map(mhc_ann.set_index('barcode').peptide)
    df_tmp['feature_type'] = 'Antigen Capture'
    df_lst.append(df_tmp)

if hsh_custom == False:
    identifier, sequence = 'Barcode','Sequence'

    df_tmp = hsh_lib.loc[hsh_lib[identifier].isin(hsh_ann.barcode), [identifier, sequence]]
    df_tmp.columns = ['id','sequence']
    df_tmp['name'] = df_tmp['id']
    df_tmp['feature_type'] = 'Antibody Capture'
    df_lst.append(df_tmp)

if mrk_custom == False:
    identifier, sequence = 'Barcode','Sequence'

    df_tmp = mrk_lib.loc[mrk_lib[identifier].isin(mrk_ann.barcode), [identifier, sequence]]
    df_tmp.columns = ['id','sequence']
    df_tmp['name'] = df_tmp['id'].map(mrk_lib.set_index(identifier).Description)
    df_tmp['feature_type'] = 'Antibody Capture'
    df_lst.append(df_tmp)

feature_df = pd.concat(df_lst, ignore_index=True)
feature_df['read'] = 'R2'
feature_df['pattern'] = '5PNNNNNNNNNN(BC)'
feature_df[['id','name','read','pattern','sequence','feature_type']].to_csv(FEATURE_REFERENCE, index=False)




