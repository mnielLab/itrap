#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
from ast import literal_eval
import seaborn as sns
import os
import yaml

import sys  
sys.path.insert(0, '../scripts')

from D_plot_specificity_matrix_utils import calc_binding_concordance

sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})

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

def plot_specificity(title, df, max_gems, save=True):
    # Sort
    df.ct = df.ct.astype(int).astype(str)
    df.sort_values(by=['epitope_rank','gems_per_specificity','binding_concordance'], #,'ct'
                       ascending=[True, False, False], inplace=True) #, True

    # devide GEMs by max concordance and outliers
    dct = df.groupby('ct').binding_concordance.max()
    df['max_conc'] = df.ct.map(dct)
    idx = df.binding_concordance == df.max_conc

    def modify_legend(h,l):
        flag = False
        labels = []
        handles = []
        for e, le in enumerate(l):
            if flag:
                #if (int(le) % 10 == 0) | (le == '1'):
                labels.append(le)
                handles.append(h[e])
            if le == 'gems_per_specificity':
                flag = True
                
        idxs = np.linspace(0,len(labels)-1,5)
        l = []
        h = []
        for i in idxs:
            l.append(labels[int(i)])
            h.append(handles[int(i)])
        return h,l #handles, labels
    
    # Style
    # https://seaborn.pydata.org/generated/seaborn.axes_style.html
    sns.set_style('ticks', {'axes.edgecolor': '0', #'axes.facecolor':'lightgrey',
                            'xtick.color': '0',
                            'ytick.color': '0'})
    sns.set_context("paper",font_scale=2)
    
    fig_height = int(df.peptide_HLA.nunique()/2) # 6
    fig = plt.figure(figsize=(20,fig_height)) # 6
    sns.scatterplot(data=df[idx], x='ct', y='peptide_HLA',
                    size='gems_per_specificity', sizes=(10,1000), size_norm=(1,max_gems),
                    hue='binding_concordance', palette='viridis_r', hue_norm=(0,1),
                    legend='full', linewidth=0)
    sns.scatterplot(data=df[~idx], x='ct', y='peptide_HLA',
                    size='gems_per_specificity', sizes=(10,1000), size_norm=(1,max_gems),
                    hue='binding_concordance', palette='viridis_r', hue_norm=(0,1),
                    legend=False, linewidth=0)
    ax = plt.gca()
    sm = plt.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=0,vmax=1), cmap='viridis_r') # hack for cbar
    sm.set_array([]) # hack for cbar
    #fig.colorbar(sm, ax=ax, pad=0.01)
    fig.colorbar(sm, ax=ax, orientation='horizontal', label='Binding Concordance', fraction=0.06*6/fig_height, pad=0.15*6/fig_height) #https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.colorbar.html

    h,l = ax.get_legend_handles_labels()
    h,l = modify_legend(h,l)
    ax.legend(h, l, bbox_to_anchor=(0.5, -0.5*6/fig_height), loc=9, frameon=False, title='GEMs', ncol=len(l))

    #plt.title('%d clonotypes (%d GEMs)' %(df.ct.nunique(), df.gem.nunique()))
    plt.xlabel('%d clonotypes (across %d GEMs)' %(df.ct.nunique(), df.gem.nunique()))
    #plt.xlabel('')
    plt.ylabel('')

    sns.despine(bottom=False, trim=True, offset={'left':-30})
    ax.set_xticks([])
    ax.set_xticklabels([])
    if save:
        plt.savefig(title, bbox_inches='tight', dpi=300)
    plt.show()


##########################################################
#                         Input                          #
##########################################################
#OPT_THR = snakemake.input.thr #'../experiments/exp13/run1_archive/cat/eval_clonotypes/threshold/opt.csv'
VALID = snakemake.input.df #'../experiments/exp13/run1_archive/cat/eval_clonotypes/valid_ct.csv'
#OS2 = '../experiments/exp13/run2/cat/eval_clonotypes/valid_ct.csv'
#HTO = '../experiments/exp13/run1_archive/brc/outfile.csv'
#GEX = snakemake.input.gex #'../experiments/exp13/run1_archive/tcr/usable_gems.txt'
#TCR = '../experiments/exp13/run1/tcr/cellranger_tot/outs/multi/vdj_t/all_contig_annotations.csv'
#TCR_ARC = '../experiments/exp13/run1_archive/tcr/cellranger_tot/outs/multi/vdj_t/all_contig_annotations.csv'


FLT = snakemake.input.flt #'indv' #'rnd_frst_oh'
IDX = snakemake.input.idx

##########################################################
#                         Output                         #
##########################################################
OUTPUT = os.path.join(os.path.dirname(snakemake.output[0]), '%s.png')

##########################################################
#                          Load                          #
##########################################################
#gex = pd.read_csv(GEX, header=None, names=['gem'])

#opt_thr = pd.read_csv(OPT_THR, index_col=0, header=None, names=['thr']).thr.dropna()

df = pd.read_csv(VALID, converters=converters)

df.fillna({'umi_count_mhc':0, 'delta_umi_mhc':0, 'umi_count_mhc_rel':0,
           'umi_count_cd8':0, 'delta_umi_cd8':0,
           'umi_count_TRA':0, 'delta_umi_TRA':0,
           'umi_count_TRB':0, 'delta_umi_TRB':0,
           'cdr3_TRA':'','cdr3_TRB':''}, inplace=True)

idx_df = pd.read_csv(IDX)

##########################################################
#                          Prep                          #
##########################################################
df = calc_binding_concordance(df.copy(), 'ct')
#df['gex'] = df.gem.isin(gex.gem)
df.single_barcode_mhc = np.where(df.single_barcode_mhc, 'pMHC singlet','pMHC multiplet')
#df['clonotype_multiplet'] = df.ct.map(df.groupby('ct').size() > 1)


##########################################################
#                         Filters                        #
##########################################################

with open(FLT, 'r') as f:
    flt = yaml.load(f, Loader=yaml.FullLoader)
globals().update(flt)
#f = open(FLT)
#globals().update(yaml.load(f)) 
#f.close()


##########################################################
#                    Compute statistics                  #
##########################################################
#df.rename(columns={'rank':'epitope_rank'}, inplace=True)

for label in labels: # do not include total in this loop
    idx = idx_df[label]
    print(label, sum(idx))
    tmp = calc_binding_concordance(df[idx].copy(), 'ct')
    
    s = '_'.join(label.split())
    
    filename = OUTPUT %s #[fn for fn in OUTPUT if label in fn]
    print(filename)
    
    max_gems = df.gems_per_specificity.max() if df.gems_per_specificity.max() < 1000 else 1000
    
    plot_specificity(filename, tmp, max_gems, save=True)
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close() # Close a figure window
    
