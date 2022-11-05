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
import seaborn as sns
from sklearn.metrics import roc_curve, auc

from D_plot_specificity_matrix_utils import (calc_binding_concordance)

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

def keep_max_conc(df):
    df = calc_binding_concordance(df.copy(), 'ct')
    assert df.binding_concordance.isna().any() == False
    
    dct = df.groupby(['ct']).binding_concordance.max()
    df['max_conc'] = df.ct.map(df.groupby(['ct']).binding_concordance.max())
    return df.binding_concordance == df.max_conc

def get_filtered_data(df, opt_thr):
    # Convert series of opt_thr to an exacutable and execute
    selection = eval(' & '.join(['(df.%s >= %d)'%(k,v) for k,v in opt_thr.iterrows()]))
    return selection

def get_multiplets(df):
    tmp = df[opt_idx & df.HLA_match_per_gem]
    dct = tmp.groupby(['ct','peptide_HLA']).gem.count() > 1
    idx = df.set_index(['ct','peptide_HLA']).index.map(dct)
    return idx.fillna(False)

def get_unique_entries(df):
    return np.where((df.cdr3_TRA==a) & (df.cdr3_TRB==b), False, True)

def get_sample_size():
    inter_entries = get_unique_entries(inter_chains)
    inter_indexes = inter_chains[inter_entries].index.to_list()
    #print(sum(inter_entries), len(inter_indexes), sum(get_unique_entries(group)))
    return min(sum(get_unique_entries(group)), len(inter_indexes))

def add_similarity_scores(ai,bi):
    mat_a = sim_tra.loc[ai, a].reset_index(drop=True).T.reset_index(drop=True).T
    mat_b = sim_trb.loc[bi, b].reset_index(drop=True).T.reset_index(drop=True).T 
    return mat_a.add(mat_b)

def avg_similarity_scores(ai,bi):
    mat_a = sim_tra.loc[ai, a].reset_index(drop=True).T.reset_index(drop=True).T.to_frame().T
    mat_b = sim_trb.loc[bi, b].reset_index(drop=True).T.reset_index(drop=True).T.to_frame().T
    return pd.concat([mat_a, mat_b], ignore_index=True).mean() #.to_frame() is important?

def get_intra_similarity(cdr3_TRAs, cdr3_TRBs, method='sum'):
    unique_entries = get_unique_entries(group)
    unique_entry_indexes = group[unique_entries].index.to_list()

    sample_size = get_sample_size() 
    # Am I cheating myself by having the chance of not sampling the most similar chains?
    # There can be quite substantial size differences: 100x
    
    sampled_idxs = random.sample(unique_entry_indexes, sample_size)
    intra_a = group.loc[sampled_idxs, 'cdr3_TRA'].values
    intra_b = group.loc[sampled_idxs, 'cdr3_TRB'].values
    
    if method == 'sum':
        combined_similarity = add_similarity_scores(intra_a, intra_b)
    else:
        combined_similarity = avg_similarity_scores(intra_a, intra_b)
    return {'score': combined_similarity.max(),
            'fraction': sum(combined_similarity > 1.8)/len(combined_similarity)}

def get_inter_similarity(cdr3_TRAs, cdr3_TRBs, method='sum'):
    # OBS! make sure the size to sample from matches the number og unique entries intra_similarity! 
    unique_entries = get_unique_entries(inter_chains)
    unique_entry_indexes = inter_chains[unique_entries].index.to_list()

    sample_size = get_sample_size() # Am I cheating myself by having the chance of not sampling the most similar chains?

    sampled_idxs = random.sample(unique_entry_indexes, sample_size)
    inter_a = inter_chains.loc[sampled_idxs, 'cdr3_TRA'].values
    inter_b = inter_chains.loc[sampled_idxs, 'cdr3_TRB'].values
    
    if method == 'sum':
        combined_similarity = add_similarity_scores(inter_a, inter_b)
    else:
        combined_similarity = avg_similarity_scores(inter_a, inter_b)
    return {'score': combined_similarity.max(),
            'fraction': sum(combined_similarity > 1.8)/len(combined_similarity)}

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

def add_number_of_observations(intra_lst, inter_lst):
    for box, lst in enumerate([intra_lst, inter_lst], start=1):
        if lst:
            median = statistics.median(lst)
            plt.text(box, median, "n: %i" %len(lst), ha='center', va='bottom')

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
    

def plot_boxplot(intra_lst, inter_lst, title='Pooled'):
    plt.figure(figsize=(3,5))
    plt.boxplot([intra_lst, inter_lst], labels=['intra', 'inter'], widths=(0.5, 0.5))
    plt.title(title)
    plt.xlim(0.6, 2.4)
    plt.ylabel("Similarity")
    
    add_number_of_observations(intra_lst, inter_lst)
    add_significance_bar((intra_lst, inter_lst))
    
    #plt.show()

    #plt.savefig(get_file('boxplot'), bbox_inches='tight')
    #plt.cla()   # Clear axis
    #plt.clf()   # Clear figure
    #plt.close()

def plot_pieplot(significant_count, total_peptides):
    plt.pie([significant_count, total_peptides-significant_count],
            labels=['significant', 'insignificant'],
            autopct=lambda p: '{:.0f} ({:.0f}%)'.format(p * total_peptides / 100, p))
    plt.title("Proportion of significant outcomes (%i)" %total_peptides)
    plt.savefig(get_file('pieplot'), bbox_inches='tight')
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close()

def add_data(scores, peptide, hue, filtering):
    n = len(scores)
    pep_lst = [peptide] * n
    hue_lst = [hue] * n
    flt_lst = [filtering] * n
    
    g = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB'])
    
    cdr3s = list(zip(g.cdr3_TRA.values, g.cdr3_TRB.values))
    cts = g.ct.values
    
    assert len(cdr3s) == n
    assert len(cts) == n
    
    tmp = pd.DataFrame(np.array([scores, pep_lst, hue_lst, flt_lst, cdr3s, cts]).T,
                       columns=['score','peptide','plateau','filtering','cdr3s','ct'])
    tmp.score = tmp.score.astype(float)
    return tmp

######################################################################################################################

# # Input
VALID = snakemake.input.valid_df
DATA = snakemake.input.original
THRESHOLD = snakemake.input.threshold
SIM_TRA = snakemake.input.sim_tra
SIM_TRB = snakemake.input.sim_trb

# # Output
#OUTPUT = snakemake.output[0]

# # Load
#valid_df = pd.read_csv(VALID, converters=converters).fillna('')
df = pd.read_csv(DATA)
opt_thr = pd.read_csv(THRESHOLD, index_col=0, header=None)
opt_idx = get_filtered_data(df, opt_thr.dropna())


#valid_ct = valid_df.ct.unqiue()

sim_tra = pd.read_csv(SIM_TRA, index_col=0).rename(index={'missing':''}, columns={'missing':''})
sim_trb = pd.read_csv(SIM_TRB, index_col=0).rename(index={'missing':''}, columns={'missing':''})
sim_tra = sim_tra[~sim_tra.index.duplicated()].copy() # Necessary?
sim_trb = sim_trb[~sim_trb.index.duplicated()].copy()

# For the similarity computation we don't want to include clonotypes that are missing a chain
df.dropna(subset=['cdr3_TRA','cdr3_TRB'], inplace=True)
# when filtering for optimal threshold it is important to have values in UMI and delta
df.fillna(value={"umi_count_mhc": 0, "delta_umi_mhc": 0, "umi_count_mhc_rel":0,
                 "umi_count_cd8": 0, "delta_umi_cd8": 0,
                 "umi_count_TRA": 0, "delta_umi_TRA": 0,
                 "umi_count_TRB": 0, "delta_umi_TRB": 0}, inplace=True) #df.fillna(0) 
df.replace('0','', inplace=True) # why?
df.sort_values(by='epitope_rank', inplace=True)
df['HLA_match_per_gem'] = df.apply(lambda row: row.HLA_mhc in row.HLA_cd8 if row.HLA_cd8 == row.HLA_cd8 else False, axis=1)

##################################################################################################
#                                          Preparations                                          #
##################################################################################################
score_dct = dict()
plt_df = pd.DataFrame(columns=['score','peptide','plateau','filtering','cdr3s','ct'])
intra_df = pd.DataFrame(columns=['score','peptide','plateau','filtering','cdr3s','ct'])
inter_df = pd.DataFrame(columns=['score','peptide','plateau','filtering','cdr3s','ct'])

filtering = ['raw','optimal threshold','HLA match', 'specificity multiplets'] #['raw','thr']
for data_type, sdf in zip(filtering, [df, df[opt_idx], df[opt_idx & df.HLA_match_per_gem], df[get_multiplets(df)]]):
    sdf = sdf[keep_max_conc(sdf)]
    for peptide, group in sdf.groupby('peptide_HLA', sort=False):
        print(peptide)
        if len(group) == 1:
            continue
        if len(group.drop_duplicates(['cdr3_TRA','cdr3_TRB'])) == 1:
            continue

        inter_chains = sdf.loc[sdf.peptide_HLA != peptide, ['cdr3_TRA', 'cdr3_TRB']]

        intra_score_peptide = list()
        inter_score_peptide = list()  

        # OBS! Maybe get the pairs directly instead of zipping them?
        # Make sure you get true pairs and not just random pairs?!
        cdr3_TRAs = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRA.values
        cdr3_TRBs = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRB.values

        assert len(cdr3_TRAs) == len(cdr3_TRBs) == len(group.loc[:,['cdr3_TRA','cdr3_TRB']].drop_duplicates())

        for index, (a,b) in enumerate(zip(cdr3_TRAs, cdr3_TRBs)):
            intra = get_intra_similarity(cdr3_TRAs, cdr3_TRBs)
            inter = get_inter_similarity(cdr3_TRAs, cdr3_TRBs)

            intra_score_peptide.append(intra['score'])
            inter_score_peptide.append(inter['score'])

        intra_df = pd.concat([intra_df, add_data(intra_score_peptide, peptide, 'intra', data_type)])
        inter_df = pd.concat([inter_df, add_data(inter_score_peptide, peptide, 'inter', data_type)])

    plt_df = pd.concat([intra_df,inter_df])
    score_dct[data_type] = list(intra_df[intra_df.filtering == data_type].score -
                                inter_df[inter_df.filtering == data_type].score)

plt_df.to_csv(snakemake.output.plt_df, index=False)

###########################
# Plot scores per peptide #
###########################
# FIX ME!
# The numbers under the plot dont match?!
#fig = plt.figure()
g = sns.catplot(x="peptide",
                y="score",
                hue="plateau",
                row="filtering",
                data=plt_df[plt_df.filtering.isin(['raw','optimal threshold'])],
                kind="box",
                aspect=2.5, sharey=True,
                legend=False)

for ax, flt in zip(g.axes, ['raw','optimal threshold']):
    add_counts(ax[0], plt_df[plt_df.filtering == flt],
               "peptide", "score",
               plt_df.peptide.unique(),
               ['intra','inter'])
    
g.set_ylabels('Similarity score')
g.set_xlabels('Peptide HLA')
g.set_xticklabels(rotation=90, ha='center')
plt.legend(bbox_to_anchor=(1.05, 1.05), loc='center', frameon=False)
sns.despine(trim=True)
plt.savefig(snakemake.output.score_per_pep, bbox_inches='tight')
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

##########################
#   Plot pooled scores   #
##########################
#fig = plt.figure()
ax = sns.boxplot(x="filtering", y="score", hue="plateau", data=plt_df)
add_significance_bar(plt_df, value='score', hue='plateau', labels='filtering')
ax.set_ylabel('Similarity score')
ax.set_xlabel('')
add_counts(ax, plt_df, "filtering", "score", filtering, ['intra','inter'])
plt.legend(bbox_to_anchor=(1.1, 0.5), loc='center', frameon=False)
sns.despine(trim=True)
plt.savefig(snakemake.output.score_pooled, bbox_inches='tight')
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

plt_df.to_csv(snakemake.output.score_pooled[:-3] + 'csv', index=False)

#############################
# Plot pooled, delta scores #
#############################
#fig = plt.figure()
comp_df = pd.DataFrame.from_dict(score_dct, orient='index').T
comp_df = pd.melt(comp_df, value_vars=['raw','optimal threshold'])
comp_df.to_csv(snakemake.output.score_pooled_delta[:-3] + 'csv', index=False)
ax = sns.boxplot(x='variable',y='value', data=comp_df)
#ax.set_xticklabels(comp_df.peptide.unique(), rotation=90, ha='center')
add_significance_bar((comp_df[comp_df.variable == 'raw'].value.dropna(),
                     comp_df[comp_df.variable == 'optimal threshold'].value.dropna()))
add_counts(ax, comp_df, "variable", "value", filtering)
ax.set_xlabel('')
ax.set_ylabel(r'$\Delta$ similarity score')
plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
plt.setp(ax.lines, color='k')
sns.despine(trim=True, bottom=False) #, offset={'bottom':20}
plt.savefig(snakemake.output.score_pooled_delta, bbox_inches='tight')
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

#############################
#            AUC            #
#############################
plt_df['label'] = np.where(plt_df.plateau == 'intra', 2, 0)

fpr = dict()
tpr = dict()
roc_auc = dict()
for i, flt in plt_df.groupby('filtering'):
    fpr[i], tpr[i], _ = roc_curve(flt.label, flt.score, pos_label=2)
    roc_auc[i] = auc(fpr[i], tpr[i])
    
plt.figure()
lw = 2
for k in fpr.keys():
    plt.plot(fpr[k], tpr[k], lw=lw, label="%s (area = %0.2f)" %(k, roc_auc[k]))
plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver operating characteristic example")
plt.legend(loc="lower right")
sns.despine(trim=True)
plt.savefig(snakemake.output.auc, bbox_inches='tight')