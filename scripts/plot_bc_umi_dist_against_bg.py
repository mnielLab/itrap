#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from ast import literal_eval
import re
from scipy import stats
from random import sample
import seaborn as sns
import matplotlib.gridspec as gridspec

#plt.style.use('ggplot')
sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})

# Input
INPUT = snakemake.input.data #'../experiments/exp13/run1/cat/tables/tcr_barcode.csv' #
#OS2 = '../experiments/exp13/run2/cat/tables/tcr_barcode.csv'

OUTPUT = snakemake.output.thrs #'../experiments/exp13/run1/cat/tables/minimum_threshold.csv' # 
PLOT = snakemake.output.plot #'umi_distributions.png' #

# Load
df = pd.read_csv(INPUT)
#os2 = pd.read_csv(OS2)

# Prep
df['umi_count_mhc'] = np.where(df.peptide_HLA.isna(), np.nan, df.umi_count_mhc)
df['umi_count_cd8'] = np.where(df.sample_id.isna(), np.nan, df.umi_count_cd8)
df['umi_count_TRA'] = np.where(df.cdr3_TRA.isna(), np.nan, df.umi_count_TRA)
df['umi_count_TRB'] = np.where(df.cdr3_TRB.isna(), np.nan, df.umi_count_TRB)

titles = {'umi_count_mhc': 'Peptide HLA barcodes',
          'umi_count_cd8': 'Hashing barcode',
          'umi_count_TRA': 'TCRa chain',
          'umi_count_TRB': 'TCRb chain'} 

threshold = dict()
for bg_factor, umis, plot in zip([df.ct.isna(), df.peptide_HLA.isna()], # & df.sample_id.isna()
                                 [('umi_count_mhc','umi_count_cd8'), ('umi_count_TRA', 'umi_count_TRB')],
                                 PLOT):
    print(sum(bg_factor))
    df['GEMs'] = np.where(bg_factor, 'Background', 'True complexes')

    # Plot
    fig, (axt, axb) = plt.subplots(2,2, figsize=(12,7), sharex=False, sharey='row')

    sns.boxplot(data=df, x=umis[0], y='GEMs', ax=axt[0])
    axt[0].set_xlabel('UMI counts')
    axt[0].set_title(titles[umis[0]])

    sns.boxplot(data=df, x=umis[1], y='GEMs', ax=axt[1])
    axt[1].set_xlabel('UMI counts')
    axt[1].set_ylabel('')
    axt[1].set_title(titles[umis[1]])

    for i, umi in enumerate(umis):
        print(umi)
        sensitivity = list()
        specificity = list()
        max_umi = int(round(df.loc[df.GEMs == 'True complexes', umi].max(), 0))
        stop = 0
        for thr in range(max_umi):
            print(thr)
            TP = len(df[(df[umi] >= thr) & (df.GEMs == 'True complexes')])
            FP = len(df[(df[umi] >= thr) & (df.GEMs == 'Background')])
            FN = len(df[(df[umi] < thr) & (df.GEMs == 'True complexes')])
            TN = len(df[(df[umi] < thr) & (df.GEMs == 'Background')])

            sensitivity.append(TP/(TP + FN))
            specificity.append(TN/(TN + FP))

            if sensitivity[-1] > specificity[-1]: # & (sensitivity[-1] > 0.95)
                print(thr, sensitivity[-1], specificity[-1])
                threshold[umi] = thr
            else: 
                if stop == 30:
                    break
                stop += 1

        df_thr = pd.DataFrame([np.arange(max_umi), sensitivity, specificity], index=['threshold','sensitivity','specificity']).T
        plt_df = pd.melt(df_thr, id_vars='threshold', value_vars=['specificity','sensitivity'])

        sns.lineplot(data=plt_df, x='threshold', y='value', hue='variable', ax=axb[i])
        axb[i].set_xlabel('Threshold')
        axb[i].set_ylabel('')

    sns.despine(trim=True, ax=axt[0])
    sns.despine(trim=True, ax=axt[1])
    sns.despine(trim=True, ax=axb[0])
    sns.despine(trim=True, ax=axb[1])
    fig.savefig(plot, bbox_inches='tight')

pd.DataFrame.from_dict(threshold, orient='index').to_csv(OUTPUT, header=False) #columns=['threshold']
