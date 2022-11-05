#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[2]:


plt.style.use('ggplot')


# # Args

# In[3]:


#EXP = "exp5"
#PLATFORM = "ILLUMINA"


# In[4]:


#MAPPING = 'KMA' # BLAST
BARCODE_SYSTEM = snakemake.params[0] #'10x' #'AKB' #


# In[5]:


if BARCODE_SYSTEM == '10x':
    BARCODE_SYSTEM_REGEX = "^(?!.*A\d+B\d+).*$"
    ANTIBODY_REGEX = "HASH"
if BARCODE_SYSTEM == 'AKB':
    BARCODE_SYSTEM_REGEX = "^A\d+B\d+"
    ANTIBODY_REGEX = "A4000"


# ## Input data

# In[6]:

map_file = snakemake.input[0]
#if MAPPING == 'KMA':
#    map_file = "/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC_" + PLATFORM + "/mapping/kma-1t1/output/mapping.clean.gz"
#if MAPPING == 'BLAST':
#    map_file = "/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC_" + PLATFORM + "/mapping/blast/blast.annotated.clean.tsv"


# OBS! We wont always have response data?!

# In[7]:


specificity_annotations = snakemake.input[1] #"/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_ILLUMINA/barcode_library/barcode_specificity_annotations.tab"
response_annotations = snakemake.input[2] #"/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp3_MHC/barcode_library/detected_responses_annotation.xlsx"


# ## Output data

# In[8]:

output = snakemake.output[0]
#if MAPPING == 'KMA':
#    output = "/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC_" + PLATFORM + "/mapping/kma-1t1/output/mapping.clean." + BARCODE_SYSTEM + ".augmented.gz"
#if MAPPING == 'BLAST':
#    output = "/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC_" + PLATFORM + "/mapping/blast/blast.annotated.clean." + BARCODE_SYSTEM + ".augmented.tsv"


# # Import input

# In[10]:


map_df = pd.read_csv(map_file) #, usecols=['query_id', 'template_id', 'gem', 'bit_score', 'alignment_length', 'tso', 'b_primer', 'anneal', 'a_primer', 'match'], sep=" ", names=["read_id", "gem", "tso", "b_primer", "anneal", "cd8_primer", "mhc_primer"]


# In[11]:


specificity_df = pd.read_excel(specificity_annotations, index_col=None, skiprows=1, names=['barcode', 'peptide', 'HLA', 'epitope'])


# In[12]:


response_df = pd.read_excel(response_annotations, index_col=None, usecols=['barcode_cd8', 'peptide'])


# ## Process mapping table

# In[13]:


map_df = map_df[(map_df.credible_alignment == True) & (map_df.barcode.str.contains(BARCODE_SYSTEM_REGEX))]


# At this point I have a table with multiple lines per GEM: each line corresponds to the best annotated read. The annotation of reads may agree on the same barcode or may disagree. Later I will count the number of reads for each barcode and only present the barcode with most reads.

# ## Annotate barcode read counts

# In[14]:


cd8_df = map_df[map_df.template_id.str.contains(ANTIBODY_REGEX, na = False)] 
mhc_df = map_df[~map_df.template_id.str.contains(ANTIBODY_REGEX, na = False)]


# In[15]:


cd8_unique_df = cd8_df.groupby(['gem'])['template_id'].apply(pd.Series.mode).to_frame().reset_index()
mhc_unique_df = mhc_df.groupby(['gem'])['template_id'].apply(pd.Series.mode).to_frame().reset_index()


# In[16]:


cd8_read_counts_df = cd8_df.groupby(['gem', 'template_id'])['query_id'].count().to_frame().reset_index()
mhc_read_counts_df = mhc_df.groupby(['gem', 'template_id'])['query_id'].count().to_frame().reset_index()
cd8_read_counts_df.rename(columns={'query_id': 'read_counts'}, inplace=True)
mhc_read_counts_df.rename(columns={'query_id': 'read_counts'}, inplace=True)


# ## Annotate specificities

# In[17]:


mhc_read_counts_df['barcode'], mhc_read_counts_df['sample'] = mhc_read_counts_df.template_id.str.rsplit("_", n=1).str


# In[18]:


mhc_read_counts_df = pd.merge(mhc_read_counts_df, specificity_df[['barcode','epitope']], how='left', on='barcode')


# ## Collapse multiple annotations per GEM into one

# In[19]:


mhc_read_counts_df.sort_values(by=['gem', 'read_counts'], inplace=True)
cd8_read_counts_df.sort_values(by=['gem', 'read_counts'], inplace=True)
mhc_read_count_diffs_df = mhc_read_counts_df.groupby(['gem']).read_counts.apply(np.array).to_frame().reset_index()
cd8_read_count_diffs_df = cd8_read_counts_df.groupby(['gem']).read_counts.apply(np.array).to_frame().reset_index()
mhc_read_count_diffs_df.rename(columns={'read_counts': 'read_counts_lst'}, inplace=True)
cd8_read_count_diffs_df.rename(columns={'read_counts': 'read_counts_lst'}, inplace=True)
mhc_read_count_diffs_df['read_count_diff'] = mhc_read_count_diffs_df.read_counts_lst.apply(lambda x: round((x[-1]-x[-2])/x[-1], 3) if len(x)>1 else 1.000)
cd8_read_count_diffs_df['read_count_diff'] = cd8_read_count_diffs_df.read_counts_lst.apply(lambda x: round((x[-1]-x[-2])/x[-1], 3) if len(x)>1 else 1.000)
mhc_read_count_diffs_df['single_barcode'] = mhc_read_count_diffs_df.read_counts_lst.apply(lambda x: True if len(x)==1 else False)
cd8_read_count_diffs_df['single_barcode'] = cd8_read_count_diffs_df.read_counts_lst.apply(lambda x: True if len(x)==1 else False)
mhc_read_count_diffs_df['template_lst'] = mhc_read_counts_df.groupby(['gem']).template_id.apply(np.array).to_frame().reset_index().template_id
cd8_read_count_diffs_df['template_lst'] = cd8_read_counts_df.groupby(['gem']).template_id.apply(np.array).to_frame().reset_index().template_id


# In[20]:


mhc_read_count_diffs_df['epitope_lst'] = mhc_read_counts_df.groupby(['gem']).epitope.apply(np.array).to_frame().reset_index().epitope


# In[21]:


cd8_mode_df = pd.merge(cd8_unique_df[['gem', 'template_id']], cd8_df, on=['gem','template_id'], how='left').merge(cd8_read_counts_df[['gem', 'template_id', 'read_counts']],
                                                                                                                  how='left', on=['gem','template_id']).merge(cd8_read_count_diffs_df, how='left', on='gem')
mhc_mode_df = pd.merge(mhc_unique_df[['gem', 'template_id']], mhc_df, on=['gem','template_id'], how='left').merge(mhc_read_counts_df[['gem', 'template_id', 'read_counts', 'epitope']],
                                                                                                                  how='left', on=['gem','template_id']).merge(mhc_read_count_diffs_df, how='left', on='gem')


# In[22]:


cd8_mode_df.sort_values(by=['gem', 'score', 'alignment_length'], inplace=True) #'credible_alignment', 'match', 
mhc_mode_df.sort_values(by=['gem', 'score', 'alignment_length'], inplace=True) #'credible_alignment', 'match', 
cd8_mode_df.drop_duplicates(subset=['gem'], keep='last', inplace=True)
mhc_mode_df.drop_duplicates(subset=['gem'], keep='last', inplace=True)


# In[23]:


print("MHC entries: %i" %mhc_mode_df.shape[0])
print("Unique GEMs: %i" %mhc_mode_df.gem.unique().shape[0])

print("CD8 entries: %i" %cd8_mode_df.shape[0])
print("Unique GEMs: %i" %cd8_mode_df.gem.unique().shape[0])



# #### OBS!
# At this point we have a table with one line per GEM. The barcode annotation is chosen from the majority vote of read counts. However, if two barcodes are annotated with equal number of reads, they are distinguished based on credible_alignment, match, bit_score, and alignment_length...
# OBS! We are not using UMI - so we don't know if the read count is inflated. And we don't take into account the alignment score or the alignment length! What if the majority of reads map poorly?

# ## Produce barcode table

# In[24]:


barcode_df = pd.merge(mhc_mode_df[['gem', 'template_id', 'template_lst', 'barcode', 'sample', 'read_counts', 'read_counts_lst', 'read_count_diff', 'single_barcode', 'credible_alignment', 'epitope_lst']],
                      cd8_mode_df[['gem', 'template_id', 'template_lst', 'barcode', 'sample', 'read_counts', 'read_counts_lst', 'read_count_diff', 'single_barcode', 'credible_alignment']],
                      how='outer', on='gem', suffixes=('_mhc', '_cd8')) #'match',


# In[25]:


assert barcode_df.shape[0] == barcode_df.gem.unique().shape[0], "Barcode dataframe was not reduced satisfyingly"


# In[26]:


print("Entries: %i" %barcode_df.shape[0])
print("Unique GEMs: %i" %barcode_df.gem.unique().shape[0])


# ## Specificity matrix

# In[27]:


specificity_matrix = mhc_mode_df.pivot(index='gem', columns='epitope', values='read_counts')


# ## Response df

# In[28]:


response_df.sort_values(by=['peptide','barcode_cd8'], inplace=True)
response_df.drop_duplicates(inplace=True)
response_df['detected_response'] = True


# ## Merge barcode and specificity

# In[29]:


barcode_specificity_df = pd.merge(barcode_df[['gem',
                                              'credible_alignment_mhc',
                                              'credible_alignment_cd8',
                                              'template_id_mhc', 'barcode_mhc', 'sample_mhc', 'read_counts_mhc', 'read_count_diff_mhc', 'single_barcode_mhc', 'read_counts_lst_mhc', 'template_lst_mhc', 'epitope_lst',
                                              'template_id_cd8', 'barcode_cd8', 'sample_cd8', 'read_counts_cd8', 'read_count_diff_cd8', 'single_barcode_cd8', 'read_counts_lst_cd8', 'template_lst_cd8']],
                                  specificity_df,
                                  how='left',
                                  left_on='barcode_mhc',
                                  right_on='barcode').merge(response_df,
                                                            how='left',
                                                            on=['barcode_cd8', 'peptide']).merge(specificity_matrix, how='left', on='gem') # 'match_mhc','match_cd8', 


# In[30]:


barcode_specificity_df['peptide_assayed'] = np.where(barcode_specificity_df.peptide.isin(response_df.peptide), True, False)


# In[31]:


print("Entries: %i" %barcode_specificity_df.shape[0])
print("Unique GEMs: %i" %barcode_specificity_df.gem.unique().shape[0])


# # Write data

# In[32]:


new_column_order = ['gem',
 'template_id_mhc',
 'umi_count_mhc',
 'single_barcode_mhc',
 'umi_count_lst_mhc',
 'read_count_lst_mhc',
 'template_lst_mhc',
 'template_id_cd8',
 'umi_count_cd8',
 'single_barcode_cd8',
 'umi_count_lst_mhc',
 'read_counts_lst_cd8',
 'template_lst_cd8',
 'detected_response',
 'peptide_assayed',
 'peptide',
 'HLA',
 'epitope',
 'epitope_lst'] + specificity_matrix.columns.to_list() #'match_mhc','match_cd8',


# In[33]:


barcode_specificity_df[new_column_order].to_csv(output, index=False, sep='\t')


# In[ ]:




