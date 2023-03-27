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
import argparse

sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})


def get_argparser():
    """Return an argparse argument parser."""
    parser = argparse.ArgumentParser(prog = 'Plot staircase',
                                     description = 'Generates staircase plots from all filters.')
    add_arguments(parser)
    return parser

def add_arguments(parser):
    parser.add_argument('--data', required=True, help='Filepath for data')
    parser.add_argument('--labels', required=True, help='Filepath for output filter labels')
    parser.add_argument('--filters', required=True, help='Filepath for filters')
    parser.add_argument('--out-dir', required=True, help='Directory to place output plots')

def HLA_cd8_converter(x):
    return x.replace("[","").replace("]","").replace(",", "").replace("'","").split(" ")

def cdr3_lst_converter(x):
    return x.replace("[","").replace("]","").replace("'","").split(" ")

def epitope_converter(x):
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
    dct = df.groupby(['ct','peptide_HLA']).gem.count() > 1
    idx = df.set_index(['ct','peptide_HLA']).index.map(dct)
    return idx.fillna(False)

def calc_binding_concordance(df, clonotype_fmt):
    gems_per_specificity        = df.groupby([clonotype_fmt,'peptide_HLA']).gem.count().to_dict()
    df['gems_per_specificity']  = df.set_index([clonotype_fmt,'peptide_HLA']).index.map(gems_per_specificity)
    gems_per_spec_hla_match     = df[df.HLA_match == True].groupby([clonotype_fmt, 'peptide_HLA']).gem.count().to_dict()
    df['gems_per_spec_hla_match'] = df.set_index([clonotype_fmt,'peptide_HLA']).index.map(gems_per_spec_hla_match)
    gems_per_clonotype          = df.groupby([clonotype_fmt]).gem.count().to_dict()
    df['gems_per_clonotype']    = df[clonotype_fmt].map(gems_per_clonotype)
    df['binding_concordance']   = df.gems_per_specificity / df.gems_per_clonotype
    df['hla_concordance']       = df.gems_per_spec_hla_match / df.gems_per_specificity
    df['hla_concordance']       = df.hla_concordance.fillna(0)
    return df

def plot_specificity(title, df, max_gems, save=True):
    # Sort
    df.ct = df.ct.astype(int).astype(str)
    df.sort_values(by=['epitope_rank','gems_per_specificity','binding_concordance'],
                       ascending=[True, False, False], inplace=True)

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
        return h,l
    
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
    sm = plt.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=0,vmax=1), cmap='viridis_r')
    sm.set_array([]) # hack for cbar
    fig.colorbar(sm, ax=ax, orientation='horizontal', label='Binding Concordance', fraction=0.06*6/fig_height, pad=0.15*6/fig_height)

    h,l = ax.get_legend_handles_labels()
    h,l = modify_legend(h,l)
    ax.legend(h, l, bbox_to_anchor=(0.5, -0.5*6/fig_height), loc=9, frameon=False, title='GEMs', ncol=len(l))

    plt.xlabel('%d clonotypes (across %d GEMs)' %(df.ct.nunique(), df.gem.nunique()))
    plt.ylabel('')

    sns.despine(bottom=False, trim=True, offset={'left':-30})
    ax.set_xticks([])
    ax.set_xticklabels([])
    if save:
        plt.savefig(title, bbox_inches='tight', dpi=100)
    plt.show()


##########################################################
#                          Main                          #
##########################################################
try:
    VALID = snakemake.input.df
    FLT = snakemake.input.lbl
    IDX = snakemake.input.flt
    OUT_DIR = os.path.dirname(snakemake.output[0])
except:
    parser = get_argparser()
    args = parser.parse_args()
    
    VALID = args.data
    FLT = args.labels
    IDX = args.filters
    OUT_DIR = args.out_dir


df = pd.read_csv(VALID, converters=converters)

df.fillna({'umi_count_mhc':0, 'delta_umi_mhc':0, 'umi_count_mhc_rel':0,
           'umi_count_cd8':0, 'delta_umi_cd8':0,
           'umi_count_TRA':0, 'delta_umi_TRA':0,
           'umi_count_TRB':0, 'delta_umi_TRB':0,
           'cdr3_TRA':'','cdr3_TRB':''}, inplace=True)

idx_df = pd.read_csv(IDX)


##########################################################
#                         Filters                        #
##########################################################

with open(FLT, 'r') as f:
    flt = yaml.load(f, Loader=yaml.FullLoader)
globals().update(flt)


##########################################################
#                    Compute statistics                  #
##########################################################
for label in labels:
    idx = idx_df[label]
    plt_df = calc_binding_concordance(df[idx].copy(), 'ct')
    
    filename = os.path.join(OUT_DIR, '%s.png' %( '_'.join(label.split()) ) )
    max_gems = df.gems_per_specificity.max() if df.gems_per_specificity.max() < 1000 else 1000
    
    plot_specificity(filename, plt_df, max_gems, save=True)
    plt.cla()
    plt.clf()
    plt.close()
