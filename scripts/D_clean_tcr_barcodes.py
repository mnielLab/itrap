#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import json

plt.style.use('ggplot')

def calc_relative_umi(grp):
    return grp.umi_count_mhc / grp.umi_count_mhc.quantile(0.9, interpolation='lower')

def calc_binding_concordance(df):
    assert df.size > 0, "df empty"
    gems_per_specificity_df = df.groupby(['clonotype','epitope']).gem.count().to_frame().reset_index()
    gems_per_specificity_df.rename(columns={'gem': 'gems_per_specificity'}, inplace=True)
    
    gems_per_clonotype_df = df.groupby(['clonotype']).gem.count().to_frame().reset_index()
    gems_per_clonotype_df.rename(columns={'gem': 'gems_per_clonotype'}, inplace=True)
    
    df = pd.merge(df, gems_per_specificity_df, on=['clonotype', 'epitope'], how='left', suffixes=('_total', '')).merge(gems_per_clonotype_df, on='clonotype', how='left', suffixes=('_total', ''))
    df['binding_concordance'] = df.gems_per_specificity / df.gems_per_clonotype
    
    return df

def check_vdjdb(credible_df):
    """
    Check if any of the clonotypes correspond to known clonotypes from VDJdb.
    The CDR3s of VDJdb may have had C removed from sequence, and thus may be a substring of the 10x CDR3 sequences.
    For alpha and beta chain separately, find the matching substring in VDJdb (if any).
    Map the VDJdb peptide annotation to our data.
    Check if the VDJ mapped peptide corresponds to the 10x annotated peptide.
    BTW, even VDJ have 'cross-reactive' TCR annotations...
    """

    # List CDR3s
    all_A3 = '|'.join(vdjdb[vdjdb.binder == 1].A3.dropna().unique())
    all_B3 = '|'.join(vdjdb[vdjdb.binder == 1].B3.dropna().unique())
    # Find matching substring in VDJdb (if any)
    credible_df['cdr3_TRA_substr'] = credible_df.cdr3_TRA.fillna('').str.findall(all_A3).apply(lambda x: x[0] if len(x)==1 else np.nan) # does empty string match anything?
    credible_df['cdr3_TRB_substr'] = credible_df.cdr3_TRB.fillna('').str.findall(all_B3).apply(lambda x: x[0] if len(x)==1 else np.nan)
    # Map
    dct = vdjdb[vdjdb.binder == 1].groupby(['A3','B3']).peptide.apply(list)
    credible_df['VDJdb_pep'] = credible_df.set_index(['cdr3_TRA_substr','cdr3_TRB_substr']).index.map(dct)
    # Check
    credible_df['VDJdb_check'] = credible_df.apply(lambda row: any([row.peptide==pep for pep in row.VDJdb_pep]) if row.VDJdb_pep==row.VDJdb_pep else np.nan, axis=1)

    return credible_df

# # Args
#PLATFORM = "IONTORRENT"
#EXP = "exp3"
#PRJ = "specificity_matrix"
#MAPPING = 'BLAST' # BLAST
#BARCODE_SYSTEM = 'AKB' #'10x'


# ## Input data
merged_annotations = snakemake.input.dat #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_CAT_IONTORRENT_KMA_AKB/tables/tcr_barcode.csv' #
HTO = snakemake.input.hto
VDJdb = snakemake.input.vdj
clone_sequencing = snakemake.params[0]
# ## Output
output = snakemake.output[0] #'/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9_CAT_IONTORRENT_KMA_AKB/tables/tcr_barcode.cleaned.csv' #
report = snakemake.output[1]

# # Import input data
hto = pd.read_csv(HTO, skiprows=1, header=None, names=['gem','seurat','umi_count_hto','feature_hto','hto_max_id','hto_sec_id',
                                                       'hto_margin','hto_classification','hto_global_class','hash_id'])
credible_df = pd.read_csv(merged_annotations) #, dtype={'single_TRA':bool, 'single_TRB':bool, 'single_chain_only':bool, 'single_tcell':bool}) # OBS! TEST THIS!
vdjdb = pd.read_csv(VDJdb)
if os.path.isfile(clone_sequencing):
    cs_df = pd.read_csv(clone_sequencing, sep='\t')
    cs_df.drop_duplicates(subset=['amino_acid'], inplace=True)
    clones = True
else:
    clones = False

rprt = {'tcr_barcode_merged_gems': list(set(credible_df.gem))}

##############################################################################################################################
#                                                            Clean                                                           #
##############################################################################################################################
# Remove NaN, None.
assert credible_df.gem.isna().sum() == 0
credible_df.dropna(subset=['gem','ct','peptide_HLA'], inplace=True) #,'template_id_mhc'
rprt['rm_nan_mhc_gems'] = list(set(credible_df.gem))
# Remove GEMs annotated with epitope 0
credible_df.drop(credible_df[credible_df.epitope == '0'].index, inplace=True)
rprt['no-0-epitope_gems'] = list(set(credible_df.gem))

###############################################################################################################################
#                                        Augment (Relative UMI count within clonotype)                                        #
###############################################################################################################################
for ct, grp in credible_df.groupby('ct'):
    umi_rel = grp.umi_count_mhc / grp.umi_count_mhc.quantile(0.9, interpolation='lower') #.max()
    credible_df.loc[umi_rel.index, 'umi_count_mhc_rel'] = umi_rel

###############################################################################################################################
#                                                     Check clones & Merge                                                    #
###############################################################################################################################
# Assign clonotypes from clone rearrangement
if clones:
    credible_df = pd.merge(credible_df, cs_df, how='left', left_on='cdr3_TRB', right_on='amino_acid')
    credible_df['clonotype'] = np.where(credible_df.sample_name.isna(), credible_df.clonotype, credible_df.sample_name)

credible_df = calc_binding_concordance(credible_df)

###############################################################################################################################
#                                                         Sanity check                                                        #
###############################################################################################################################
# Only include GEMs where there is a single clonotype 
#if not credible_df.groupby(['gem']).clonotype.nunique().eq(1).all():
#    print(credible_df.groupby(['gem']).clonotype.nunique().eq(1))
if not credible_df.groupby(['gem']).ct.nunique().eq(1).all():
    print(credible_df.groupby(['gem']).ct.nunique().eq(1))

print("Cleaned df")
print("Rows: %i" %credible_df.shape[0])
print("GEMs: %i" %credible_df.gem.unique().shape[0])

###############################################################################################################################
#                                                         Check VDJdb                                                         #
###############################################################################################################################
credible_df = check_vdjdb(credible_df)

###############################################################################################################################
#                                                       ADD HTO analysis                                                      #
###############################################################################################################################
df = pd.merge(credible_df, hto, how='left', on='gem')

###############################################################################################################################
#                                                            Write                                                            #
###############################################################################################################################
#credible_df.to_csv(output, index=False)
df.to_csv(output, index=False)

with open(report, 'w') as outfile:
    json.dump(rprt, outfile)


# # Extra?
# Should I only include GEMs with both a TRA and a TRB?
# 
# Should I exclude GEMs where a TCR chain could not be annotated unambiguously?

# In[ ]:




