#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

#plt.style.use('ggplot')
sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})
sns.set_context("paper",font_scale=2)
mpl.rcParams['axes.grid'] = False

# # Load
df = pd.read_csv(snakemake.input[0])
df['pMHC category'] = np.where(df.multiplets_mhc, 'multiplet','singlet')

# Extract minimum UMI count from output filename
label = snakemake.params.sorting

if 'tot' not in label:
    l = 1 if 'pos' in label else 0
    df = df[df.label == l].copy()
    
#min_umi = int(re.search('b(\d+).', snakemake.output[0]).group(1))
#print(min_umi)
#
#df = credible_df[credible_df.umi_count_mhc >= min_umi]

matrix = (df.fillna(0).groupby(['pMHC category','tcr_category']).gem.size().to_frame().reset_index()
          .pivot(index='pMHC category', columns='tcr_category', values='gem'))

fig, ax = plt.subplots() #figsize=(2.4,3.6) #figsize=(10,10)
im = ax.imshow(matrix, interpolation='nearest', cmap='viridis')

# Loop over data dimensions and create text annotations.
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        text = ax.text(j, i, "%.f" %matrix.iloc[i, j],
                       ha="center", va="center", color="w")

ax.set_xticks(np.arange(matrix.shape[1]))
ax.set_yticks(np.arange(matrix.shape[0]))
       
ax.set_xticklabels(['\n'.join(l.split()) for l in matrix.columns]) #, rotation=40, ha='right' #sorted(df.tcr_category.value_counts().index)
ax.set_yticklabels(matrix.index, rotation=90, va='center') #sorted(df.multiplets_mhc.value_counts().index)
#
# HACK: with matplotlib==3.1.1 the yaxis is cropped..?
ax.set_ylim(sorted((-0.5, matrix.shape[0]-0.5), reverse=True))

plt.xlabel('TCR chain annotations') #A- and B-c
plt.ylabel('pMHC annotations')

sns.despine(left=True,top=True,bottom=True,right=True)

#plt.title('TCR chain annotations & BC multiplets\n(in %i GEMs with min. %i BC UMIs)' %(matrix.sum().sum(), min_umi))
plt.savefig(snakemake.output[0], bbox_inches='tight', dpi=300)