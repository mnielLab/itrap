#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from ast import literal_eval
import re
from scipy import stats
from random import sample
#import scipy.stats as stats
import random
import re
import statistics
from scipy.stats import wilcoxon
from sklearn.metrics import roc_curve, auc

#plt.style.use('ggplot')

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


import sys  
sys.path.insert(0, '../scripts')

#from D_plot_specificity_matrix_utils import (peptide_per_clonotype_by_gem_size,
#                                             multiple_peptides_per_gem_w_filtering,
#                                             calc_binding_concordance,
#                                             epitope_sorter_index,
#                                             peptides_per_gem)

from D_plot_specificity_matrix_utils import calc_binding_concordance

import seaborn as sns
#sns.set_style('white')
#sns.set_context('paper') #, font_scale=2
sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})

###########################################################
#                        Functions                        #
###########################################################
def paired_t_test(x1 ,x2):
    assert len(x1) == len(x2)
    statistic, pvalue = stats.ttest_rel(x1, x2)
    if (pvalue/2.0 < 0.05) & (statistic > 0) & (len(x1) > 9):
        return {'test':True, 'pvalue':pvalue}
    else:
        return {'test':False, 'pvalue':pvalue}

def t_test(x1, x2):
    statistic, pvalue = stats.ttest_ind(x1, x2, equal_var=False, nan_policy='omit', alternative='less')
    if (pvalue < 0.05) & (statistic < 0) & (len(x1) > 9):
        return {'test':True, 'pvalue':pvalue}
    else:
        return {'test':False, 'pvalue':pvalue}
    
def wilcoxon_test(intra, inter):
    w, p = wilcoxon(np.array(intra)-np.array(inter), alternative='greater')
    if (p < 0.05) & (len(intra) > 9):
        return {'test':True, 'pvalue':p}
    else:
        return {'test':False, 'pvalue':p}

def add_significance_bar(data, hue=False, value=None, labels=None):
    print('adding significance bars')
    h1 = 1.02
    h2 = 1.025
    h3 = 1.03
    h4 = 1.035
    
    if hue:
        print('with hue')
        hue1, hue2 = data[hue].unique()
        
        y0 = data[value].max()
        y1 = y0 * h1
        y2 = y0 * h2
        y3 = y0 * h3
        y4 = y0 * h4
        
        for i,l in enumerate(data[labels].unique()):
            intra_lst = data.loc[(data[labels] == l) & (data[hue] == hue1), value]
            inter_lst = data.loc[(data[labels] == l) & (data[hue] == hue2), value]
            
            t = paired_t_test(intra_lst, inter_lst)
            print(t)
            if t['test'] and t['pvalue'] < 0.05:
                pass
            else:
                return
            
            x1, x2 = i-0.2, i+0.2
            plt.plot([x1,x1,x2,x2], [y1,y2,y2,y1], lw=1.5, c='k')
            plt.text(i, y3, "p = %.2e" %t['pvalue'], ha='center', va='bottom', color='k')
            plt.plot(1, y4)
            
    else:
        print('without hue')
        intra_lst, inter_lst = data
        if len(intra_lst) == len(inter_lst):
            t = paired_t_test(intra_lst, inter_lst)#['pvalue']
            x = [1,1,2,2]
        else:
            t = t_test(intra_lst, inter_lst)
            x = [0,0,1,1]
        print(t)
        if t['test'] and t['pvalue'] < 0.05:
            pass
        else:
            return

        y0 =  max(max(intra_lst), max(inter_lst))
        y1 = y0 * h1
        y2 = y0 * h2
        y3 = y0 * h3
        y4 = y0 * h4
        
        print('plotting')
        plt.plot(x, [y1,y2,y2,y1], lw=1.5, c='k')
        plt.text(np.mean(x), y3, "p = %.2e" %t['pvalue'], ha='center', va='bottom', color='k')
        plt.plot(1, y4)
    
def add_counts(ax, plt_df, x_col, y_col, order, hue=None):
    print('add counts')
    if hue is None:
        d = 1
    else:
        d = len(hue)
    
    y_lim = ax.get_ylim()
    y_rng = (y_lim[1] - y_lim[0]) * 0.01
    ax.set_ylim(y_lim[0]-2*y_rng, y_lim[1])
    y_min = round(ax.get_ylim()[0], 2) #y_lim[0]-y_rng
    x_pos = np.arange(len(order))
    y_pos = pd.Series([y_min]*len(order), index=order)
    counts = plt_df.dropna()[x_col].value_counts()/d
    counts = counts.reindex(order)
    
    if len(order) > 5:
        txt = '{:.0f}'
    else:
        txt = 'N={:.0f}'

    for p,n,m in zip(x_pos,counts,y_pos):
        print(p,n,m)
        if not np.isnan(m):
            ax.annotate(txt.format(n), xy=(p, m), xycoords='data', ha='center', va='bottom')

##########################################################
#                         Inputs                         #
##########################################################
SIM_INPUT = snakemake.input.sim
AUC_INPUT = snakemake.input.auc

filter_set = snakemake.params.flt
##########################################################
#                         Output                         #
##########################################################


##########################################################
#                          Load                          #
##########################################################
test_df = pd.read_csv(SIM_INPUT).dropna()
auc_df = pd.read_csv(AUC_INPUT)

labels = test_df.filtering.unique()
print(labels)
palette = auc_df.set_index('filtering').palette.to_dict()

# Computing significance of separation between intra and inter for each random sample
sign_add = pd.Series([0]*len(labels), index=labels)
for s in test_df.rnd_sample.unique():
    lst_scores = test_df[test_df.rnd_sample == s].groupby(['filtering','plateau'],sort=False).score.apply(list).to_frame().reset_index() #
    #lst_scores = test_df.dropna().groupby(['filtering','plateau'],sort=False).score.apply(list).to_frame().reset_index()
    lst_lst_scores = lst_scores.groupby('filtering', sort=False).score.apply(list)
    # by not sorting during grouping I ensure that intra scores comes first (x[0]) and inter comes last (x[-1])
    sign_test = lst_lst_scores.apply(lambda x: wilcoxon_test(x[0],x[-1])['test']) #wilcoxon_test()
    sign_add = sign_add + sign_test
# If more than half of the tests were significant I'll mark it as significant
sign_test = sign_add > len(test_df.rnd_sample.unique())/2

##########################################################
#                 Plot similarity metrics                #
##########################################################
sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})
sns.set_context("paper",font_scale=1.5)


###########################
# Plot scores per peptide #
###########################
selected_filters = ['total','optimal threshold']
g = sns.catplot(y="peptide",
                x="score",
                hue="plateau", palette=['white','grey'],
                col="filtering",
                data=test_df[test_df.filtering.isin(selected_filters)],
                kind="box",sharey=True,
                legend=False, linewidth=0.5, fliersize=1, aspect=1.2)

for ax, flt in zip(g.axes[0], selected_filters):
    ct_per_pep = test_df[(test_df.rnd_sample==0) &
                         (test_df.plateau == 'intra') & 
                         (test_df.filtering == flt)].groupby(['peptide'],sort=False).size() #.dropna()
    
    for i,pep in enumerate(test_df.peptide.unique()):
        if pep in ct_per_pep.index:
            n = ct_per_pep[pep]
            ax.annotate(f'{n}', xy=(ax.get_xlim()[0],i), xytext=(-10, -1), textcoords='offset points', va='center', ha='right')

g.set_xlabels('Similarity score')
g.set_ylabels('') #Peptide HLA
#g.set_xticklabels(rotation=90, ha='center')
#plt.xticks(rotation=90, ha='center')
plt.legend(bbox_to_anchor=(0.5, -0.1), loc='center', frameon=False, ncol=2, bbox_transform=g.fig.transFigure)
#plt.legend(bbox_to_anchor=(1.2, 0.5), loc=6, frameon=False)
sns.despine(trim=True, offset={'left':55})
plt.subplots_adjust(wspace=0.40) #38

plt.savefig(snakemake.output.score_per_pep, bbox_inches='tight', dpi=300)
#plt.savefig('../experiments/exp13/run1_archive/plt/eval_filters/sim.%s.png' %filter_set, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window


##########################
#   Plot pooled scores   #
##########################
#fig = plt.figure()
ax = sns.boxplot(y="filtering", x="score", hue="plateau", palette=['white','grey'], data=test_df, linewidth=0.5, fliersize=1)
#add_significance_bar(plt_df, value='score', hue='plateau', labels='filtering')
ax.set_xlabel('Similarity score')
ax.set_ylabel('')

# Annotate counts
ct_per_flt = test_df[(test_df.rnd_sample==0) &
                     (test_df.plateau == 'intra')].dropna().groupby(['filtering'],sort=False).size()
for i,n in enumerate(ct_per_flt):
    ax.annotate(f'{n}', xy=(ax.get_xlim()[0],i), xytext=(15, -1), textcoords='offset points', va='center', ha='right')

# Add significance bar
for i, c in enumerate(labels):
    if sign_test[c]:
        y = [i+0.22, i+0.22, i-0.22, i-0.22] #[0,0,1,1] #
        x0 = 2
        x1 = x0 * 1.007
        x2 = x0 * 1.010
        x3 = x0 * 1.020

        ax.plot([x1, x2, x2, x1], y, lw=0.7, c='0') #lw=1.5, 
        ax.plot(x3, np.mean(y), marker="*", c='0')
            
#for i, c in enumerate(labels):
#    ann = '*' if sign_test[c] else ''
#    ax.annotate(ann, xy=(2, i), #p.get_width()
#        xytext=(5, -1), textcoords='offset points', ha="left", va="center", size=12)

plt.legend(bbox_to_anchor=(0.5, -0.3), ncol=2, loc=10, frameon=False)
sns.despine(trim=True, offset={'left':30})
plt.savefig(snakemake.output.score_pooled, bbox_inches='tight', dpi=300)
#plt.savefig('../experiments/exp13/run1_archive/plt/eval_filters/sim.pool.%s.png' %filter_set, bbox_inches='tight',dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

##########################
#        Plot AUC        #
##########################
sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})
sns.set_context("paper",font_scale=2)

sns.barplot(data=auc_df, hue='filtering',y='value',x='variable', palette=palette, ci='sd')
plt.legend(bbox_to_anchor=(1.01, 0.5), loc=6, frameon=False)
plt.ylabel('')
plt.xlabel('')
sns.despine()
#plt.xticks(rotation=90)
plt.savefig(snakemake.output.auc, bbox_inches='tight', dpi=300)
#plt.savefig('../experiments/exp13/run1_archive/plt/eval_filters/auc.bar.%s.png' %filter_set, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

##########################
#        Plot ROC        #
##########################
fpr = dict()
tpr = dict()
roc_auc = dict()
for i, flt in test_df.dropna().groupby('filtering'):
    fpr[i], tpr[i], _ = roc_curve(flt.label, flt.score, pos_label=2)
    roc_auc[i] = auc(fpr[i], tpr[i])

plt.figure()
lw = 2
for k,c in zip(labels,palette):
    plt.plot(fpr[k], tpr[k], c=palette[k], lw=lw, label="%s (area = %0.2f)" %(k, roc_auc[k]))
plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.gca().set_aspect('equal', 'box')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
#plt.title("Receiver operating characteristic example")
plt.legend(bbox_to_anchor=(1.1, 0.5), loc=6, frameon=False)
sns.despine(trim=True)
plt.savefig(snakemake.output.roc, bbox_inches='tight', dpi=300)
#plt.savefig('../experiments/exp13/run1_archive/plt/eval_filters/auc.roc.%s.png' %filter_set, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window


