#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
from ast import literal_eval
import seaborn as sns
import yaml
import decimal

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


##########################################################
#                         Input                          #
##########################################################
VALID = snakemake.input.df #'../experiments/exp13/run1_archive/cat/eval_clonotypes/valid_ct.csv'
FLT = snakemake.input.flt #'indv' #'rnd_frst_oh'
IDX = snakemake.input.idx

AUC = snakemake.input.auc #'tmp_files/similarity.auc.%s.csv' %filter_set


##########################################################
#                          Load                          #
##########################################################
auc_df = pd.read_csv(AUC)

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

labels = labels[1:]
##########################################################
#                          Prep                          #
##########################################################
df = calc_binding_concordance(df.copy(), 'ct')
#df['gex'] = df.gem.isin(gex.gem)
df.single_barcode_mhc = np.where(df.single_barcode_mhc, 'pMHC singlet','pMHC multiplet')
df['clonotype_multiplet'] = df.ct.map(df.groupby('ct').size() > 1)

total_gems = len(df.gem.unique())
total_cts = len(df.ct.unique())
#total_acc = round(df.train_label.sum() / len(df.train_label.dropna()) * 100, 1)
total_acc = round(df.pep_match.sum() / len(df.pep_match.dropna()) * 100, 1)
total_conc = round(calc_binding_concordance(df[df.clonotype_multiplet].copy(), 'ct').binding_concordance.mean() * 100, 1)



##########################################################
#                       Prep plot df                     #
##########################################################

plt_df = pd.DataFrame(columns=['counts','percent','variable','filters'])
plt_df.loc[0,['counts','percent','variable','filters']] = total_gems, 100, 'GEMs', 'total'
plt_df.loc[1,['counts','percent','variable','filters']] = total_cts, 100, 'clonotypes', 'total'
plt_df.loc[2,['counts','percent','variable','filters']] = total_acc, total_acc, 'accuracy', 'total'
plt_df.loc[3,['counts','percent','variable','filters']] = total_conc, total_conc, 'avg. conc.', 'total'

# ### Initializing dataframe for TCR category summary
data = ['total'] + labels 
binding_category = ['normal','outlier']
tcr_category = ['missing chain','multiple chains','unique chains']

indexes = len(binding_category)*len(tcr_category)*len(data)

# Initiate dataframe
smry = pd.DataFrame(columns=['freq','frac','pos'],index=np.arange(indexes))
smry.index = pd.MultiIndex.from_product([data, binding_category, tcr_category],
                                        names=['data', 'binding_category','tcr_category'])


# #### Populating dataframe with summary of total data
tmp = calc_binding_concordance(df.copy(), 'ct')
tmp['binding_category'] = np.where(tmp.binding_concordance > 0.5, 'normal','outlier')

stmp = tmp.groupby(['binding_category','tcr_category']).size().to_frame().rename(columns={0:'freq'})
stmp.loc['normal','frac'] = stmp.loc['normal','freq'].values / stmp.loc['normal','freq'].sum() * 100
stmp.loc['outlier','frac'] = stmp.loc['outlier','freq'].values / stmp.loc['outlier','freq'].sum() * 100
stmp['pos'] = 0

stmp = pd.concat({'total': stmp}, names=['data']) #stmp['data'] = 'total'

smry.loc[stmp.index, ['freq','frac','pos']] = stmp


# ### Initializing dataframe for pMHC summary
data = ['total'] + labels
binding_category = ['normal','outlier']
pmhc_category = ['pMHC multiplet','pMHC singlet']

indexes = len(binding_category)*len(pmhc_category)*len(data)

# Initiate dataframe
pmhc = pd.DataFrame(columns=['freq','frac','pos'], index=np.arange(indexes))
pmhc.index = pd.MultiIndex.from_product([data, binding_category, pmhc_category],
                                        names=['data', 'binding_category','single_barcode_mhc'])


# #### Populating dataframe with summary of total data
ptmp = tmp.groupby(['binding_category','single_barcode_mhc']).size().to_frame().rename(columns={0:'freq'})

ptmp.loc['normal','frac'] = ptmp.loc['normal','freq'].values / ptmp.loc['normal','freq'].sum() * 100
ptmp.loc['outlier','frac'] = ptmp.loc['outlier','freq'].values / ptmp.loc['outlier','freq'].sum() * 100
ptmp = pd.concat({'total': ptmp}, names=['data']) #ptmp['data'] = 'total'

pmhc.loc[ptmp.index, ['freq','frac']] = ptmp

##########################################################
#                    Compute statistics                  #
##########################################################
df.rename(columns={'rank':'epitope_rank'}, inplace=True)

i = 4

#for idx, label in zip(filterings,labels): # do not include total in this loop
for label in labels:
    idx = idx_df[label]
    tmp = calc_binding_concordance(df[idx].copy(), 'ct')
    tmp['binding_category'] = np.where(tmp.binding_concordance > 0.5, 'normal','outlier')
    
    n_g = len(df[idx].gem.unique())
    n_c = len(df[idx].ct.unique())
    #n_a = round(df[idx].train_label.sum() / len(df[idx].train_label.dropna()) * 100, 1)
    n_a = round(df[idx].pep_match.sum() / len(df[idx].pep_match.dropna()) * 100, 1)
    n_s = round(tmp[tmp.clonotype_multiplet].binding_concordance.mean() * 100, 1)
    plt_df.loc[i , ['counts','percent','variable','filters']] = n_g, n_g / total_gems * 100, 'GEMs', label
    plt_df.loc[i+1,['counts','percent','variable','filters']] = n_c, n_c / total_cts * 100, 'clonotypes', label
    plt_df.loc[i+2,['counts','percent','variable','filters']] = n_a, n_a, 'accuracy', label
    plt_df.loc[i+3,['counts','percent','variable','filters']] = n_s, n_s, 'avg. conc.', label
    
    stmp = tmp.groupby(['binding_category','tcr_category']).size().to_frame().rename(columns={0:'freq'})
    #smry = smry.reset_index().groupby('tcr_category').freq.sum().to_frame()
    stmp.loc['normal','frac'] = stmp.loc['normal','freq'].values / stmp.loc['normal','freq'].sum() * 100
    stmp.loc['outlier','frac'] = stmp.loc['outlier','freq'].values / stmp.loc['outlier','freq'].sum() * 100
    stmp['pos'] = i / 4
    stmp = pd.concat({label: stmp}, names=['data']) #stmp['data'] = label
    # Fill
    smry.loc[stmp.index, ['freq','frac','pos']] = stmp
    
    ptmp = tmp.groupby(['binding_category','single_barcode_mhc']).size().to_frame().rename(columns={0:'freq'})
    ptmp.loc['normal','frac'] = ptmp.loc['normal','freq'].values / ptmp.loc['normal','freq'].sum() * 100
    ptmp.loc['outlier','frac'] = ptmp.loc['outlier','freq'].values / ptmp.loc['outlier','freq'].sum() * 100
    ptmp = pd.concat({label: ptmp}, names=['data']) #ptmp['data'] = label
    pmhc.loc[ptmp.index, ['freq','frac']] = ptmp
    
    i += 4
    

auc_df.rename(columns={'filtering':'filters'}, inplace=True)
auc_df['percent'] = auc_df.value * 100
dct = auc_df.groupby(['filters','variable']).percent.mean()
auc_df['counts'] = auc_df.set_index(['filters','variable']).index.map(dct)


##########################################################
#                           AUC                          #
##########################################################
auc_df.filters = np.where(auc_df.filters == 'is functional cell', 'is viable cell', auc_df.filters)
test_df = pd.concat([plt_df, auc_df])
test_df = test_df[test_df.variable.isin(['GEMs','accuracy','avg. conc.','AUC','AUC 0.1'])]
test_df['tmp'] = np.where(test_df.variable == 'AUC', test_df['counts'], test_df['percent'])

tmp = test_df.drop_duplicates(subset=['variable','filters']).groupby(['filters'], sort=False).tmp.mean().to_frame().reset_index()
tmp['variable'] = 'summary'
tmp['counts'] = tmp.tmp
tmp['percent'] = tmp.tmp

test_df = pd.concat([test_df, tmp])
test_df = test_df[test_df.variable.isin(['GEMs','accuracy','avg. conc.','AUC','summary'])].copy()
test_df.variable = np.where(test_df.variable == 'clonotypes','clones',test_df.variable)
test_df.variable = np.where(test_df.variable == 'avg. conc.','avg.\nconc.',test_df.variable)

cnt = test_df.groupby('filters', sort=False).counts.unique()
cnt_flat = [x for y in cnt.to_list() for x in y]

sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})
sns.set_context("paper",font_scale=2)

##########################################################
#            Plot statistics (bar) + AUC + summary       #
##########################################################
plt.figure(figsize=(10,4))#7
ax = sns.barplot(data=test_df, hue='filters',y='percent',x='variable', ci='sd',
                 palette=palette) #['grey','#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494']
for p, c in zip(ax.patches, cnt_flat):
    yoff = 7
    if type(c) == int:
        ax.annotate("%d" % c, xy=(p.get_x()+p.get_width()/2, p.get_height()), #p.get_width()
            xytext=(2, yoff), textcoords='offset points', ha="center", va="bottom", size=8, rotation=90)
    else:
        if decimal.Decimal(str(c)).as_tuple().exponent < -2:
            yoff = 12
        ax.annotate("%.1f" % round(c,1), xy=(p.get_x()+p.get_width()/2, p.get_height()), #p.get_width()
                xytext=(2, yoff), textcoords='offset points', ha="center", va="bottom", size=8, rotation=90) #
    
ax.legend(title=False, frameon=False, bbox_to_anchor=(0.5,-0.2), loc='upper center')
#ax.legend(title=False, frameon=False, bbox_to_anchor=(1.01,0.5), loc='center left', prop={'size': 9})
sns.despine()#trim=True, offset={'left':10}
plt.xlabel('')
plt.ylabel('')
plt.savefig(snakemake.output.sco, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window
#plt.savefig('../experiments/exp13/run1_archive/plt/eval_filters/sco.%s.png'%filter_set, bbox_inches='tight',dpi=300)
#plt.show()




##########################################################
#                   Plot statistics (bar)                #
##########################################################
smry.fillna(0, inplace=True)
smry.reset_index(inplace=True)
#smry['data'] = smry.index.get_level_values('data')
pmhc.fillna(0, inplace=True)
pmhc.reset_index(inplace=True)

cnt = plt_df.groupby('filters', sort=False).counts.apply(list) #'filters'
prc = plt_df.groupby('filters', sort=False).percent.apply(list)
cnt_flat = [x for y in cnt.to_list() for x in y] #[item for sublist in t for item in sublist]

plt.figure(figsize=(13,5))#7
ax = sns.barplot(data=plt_df, hue='filters',y='percent',x='variable',
                 palette=palette) #['grey','#ffffcc','#c7e9b4','#7fcdbb','#41b6c4','#2c7fb8','#253494']
for p, c in zip(ax.patches, cnt_flat):
    if type(c) == int:
        ax.annotate("%d" % c, xy=(p.get_x()+p.get_width()/2, p.get_height()), #p.get_width()
            xytext=(0, 5), textcoords='offset points', ha="center", va="bottom", size=8, rotation=90)
    else:
        ax.annotate("%.1f" % c, xy=(p.get_x()+p.get_width()/2, p.get_height()), #p.get_width()
                xytext=(0, 5), textcoords='offset points', ha="center", va="bottom", size=8, rotation=90) #
    
ax.legend(title=False, frameon=False, bbox_to_anchor=(1.05,0.5), loc=6)
sns.despine()#trim=True, offset={'left':10}
plt.xlabel('')
plt.ylabel('')
#plt.savefig('../experiments/exp13/run1_archive/plt/eval_filters/sco.indv_thr.png', bbox_inches='tight')
#plt.show()
plt.savefig(snakemake.output.old, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window


##########################################################
#                   TCR category ouliters                #
##########################################################
plt.figure(figsize=(7,5))
g = sns.catplot(data=smry.reset_index(), col='binding_category', y='frac', x='tcr_category', hue='data',
                kind='bar', sharey=True, margin_titles=True, palette=palette, legend=False)

for ax, d in zip(g.axes.flatten(), smry.binding_category.unique()): #smry.index.get_level_values(1).unique()
    for p, c in zip(ax.patches, smry[smry.binding_category == d].freq.values): #.sort_values(by=['pos','tcr_category'])
        if any(smry.freq > 1000):
            ax.annotate("%d" % c, xy=(p.get_x()+p.get_width()/2, p.get_height()), #
                        xytext=(1, 5), textcoords='offset points', ha="center", va="bottom", size=8, rotation=90)
        else:
            ax.annotate("%d" % c, xy=(p.get_x()+p.get_width()/2, p.get_height()), #p.get_width()
                xytext=(0, 5), textcoords='offset points', ha="center", va="center", size=8)

g.set_titles(col_template="{col_name}")#, row_template="{row_name}"
g.set_xlabels('')
g.set_ylabels('')
#plt.suptitle('Normalized difference between outliers and "normal" observations', y=1.05)
sns.despine() #trim=True, offset={'left':10}
#plt.savefig('../experiments/exp13/run1_archive/plt/eval_filters/tcr.thr_vs_hto_gex.png')
#plt.show()
plt.savefig(snakemake.output.tcr, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window


##########################################################
#                   pMHC category ouliters               #
##########################################################
plt.figure(figsize=(5,7))
g = sns.catplot(data=pmhc.reset_index(), col='binding_category', y='frac', x='single_barcode_mhc',
                hue='data', kind='bar', sharey=True, margin_titles=True, palette=palette, legend=False)

for ax, d in zip(g.axes.flatten(), pmhc.binding_category.unique()):
    for p, c in zip(ax.patches, pmhc[pmhc.binding_category == d].freq.values):
        if any(pmhc.sort_index().freq.values > 1000):
            ax.annotate("%d" % c, xy=(p.get_x()+p.get_width()/2, p.get_height()), #p.get_width()
                xytext=(1, 5), textcoords='offset points', ha="center", va="bottom", size=10, rotation=90)
        else:
            ax.annotate("%d" % c, xy=(p.get_x()+p.get_width()/2, p.get_height()), #p.get_width()
                xytext=(0, 5), textcoords='offset points', ha="center", va="center", size=10)
        
g.set_titles(col_template="{col_name}")#, row_template="{row_name}"
g.set_xlabels('')
g.set_ylabels('')
#plt.suptitle('Single barcode pMHC annotations', y=1.05)
#plt.savefig('../experiments/exp13/run1_archive/plt/eval_filters/pep.thr_vs_hto_gex.png')
#plt.show()
plt.savefig(snakemake.output.pep, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window
