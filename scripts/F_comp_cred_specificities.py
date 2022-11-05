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

import os
import sys
from D_plot_specificity_matrix_utils import (calc_binding_concordance)

sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})
sns.set_context("paper",font_scale=2)

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
              'umi_count_lst_cd8': literal_converter,
              'umi_count_lst_TRA': literal_converter,'umi_count_lst_TRB': literal_converter,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'cdr3_lst_TRB': cdr3_lst_converter,
              'HLA_lst_mhc': cdr3_lst_converter,
              'HLA_pool_cd8':cdr3_lst_converter,
              'HLA_cd8': HLA_cd8_converter,
              'HLA_lst_cd8':literal_converter,'sample_id_lst':literal_converter} #

#########################################################################################################
#                                                 Class                                                 #
#########################################################################################################
class Evaluate_Clonotype():
    """
    An instance is a clonotype subset of the data.
    """
    value_bin = set() # clonotypes Evaluate_Clonotype.value_bin
    trash_bin = set() # peptide-HLAs
    ct_checks = dict()
    
    def __init__(self, df, ct, selected_clonotypes, use_relative_umi=False, variable='peptide_HLA'):
        self.df = df[df.ct == ct]
        self.ct = int(ct)
        self.idx = self.df.index
        self.rel_umi = use_relative_umi
        self.fig_flg = 'None'
        self.variable = variable
        
        if self.variable == 'peptide_HLA':
            self.var_lst = 'peptide_HLA_lst'
            self.umi_lst = 'umi_count_lst_mhc'
        elif self.variable == 'sample_id':
            self.var_lst = 'sample_id_lst'
            self.umi_lst = 'umi_count_lst_cd8'
        elif self.variable == 'HLA_cd8':
            self.var_lst = 'HLA_pool_cd8'
            self.umi_lst = 'umi_count_lst_cd8'
        
        # Initialize list of clonotypes
        self.sel_cts = selected_clonotypes
        self.sel_cts.index = self.sel_cts.index.astype(int)
        
        # Initialize matrix of query count per GEM (not method (df is not self <- WHY?))
        self.queries = df[self.var_lst].explode().drop_duplicates() #.unique()
        self.mat = pd.DataFrame(index=self.queries, columns=df.gem.unique())
        
        # Count no. of GEMs within grp that are annotated with a specific pMHC (method?)
        self.gems_per_query = self.df.explode(self.variable).groupby(self.variable).size() # explode only really necessary for HLA_pool_cd8, but doesnt make a difference for other variables.
        self.gems_per_all_q = pd.concat([self.gems_per_query,
                                         pd.Series(0, index=self.queries[~self.queries.isin(self.gems_per_query.index)])])
        
    def sum_umi(self):
        # Sum UMIs for each peptide across GEMs (multiple peptides per GEM)
        for idx, row in self.df.iterrows():
            if self.variable == 'HLA_cd8':
                var = 'HLA_lst_cd8'
                var_lst = [item for sublist in row[var] for item in sublist if item != '']
                umi_lst = [row[self.umi_lst][i] for i, sublist in enumerate(row[var]) for item in sublist if item != '']
            else:
                var = self.var_lst
                var_lst = row[self.var_lst]
                umi_lst = row[self.umi_lst]
                
            if len(row[self.umi_lst]) == len(row[var]):
                self.mat.loc[var_lst, row.gem] = umi_lst
            else:
                self.mat.loc[var_lst, row.gem] = [0] * len(var_lst)
        self.mat.fillna(0, inplace=True) # outcomented..
                
    def calc_summary(self):
        self.summary_df = self.mat.sum(axis=1).sort_values(ascending=False).to_frame().rename(columns={0:'s'})
        self.summary_df['avg'] = self.mat.mean(axis=1)
        self.summary_df['col'] = 'grey'
        self.summary_df['r'] = self.summary_df.s / self.summary_df.s.max() # Unnecessary
        
    def calc_relative_umi(self):
        if self.variable == 'peptide_HLA':
            umi = 'umi_count_mhc'
        else:
            umi = 'umi_count_cd8'
        if self.rel_umi:
            self.mat = self.mat / self.df[umi].quantile(0.9, interpolation='lower')
        return self.df[umi] / self.df[umi].quantile(0.9, interpolation='lower')
    
    #def select_queries(self, n=11):
    #    self.selected_queries = self.summary_df.head(n).index
    #    self.selected_mat = self.mat[self.mat.index.isin(self.selected_queries)]
        
    def select_queries(self, n=11):
        self.selected_queries = self.summary_df.head(n).index
        self.selected_mat = self.mat.loc[self.mat.index.isin(self.selected_queries), self.df.gem]

    #def transform_data_for_plotting(self):
    #    """
    #    For each row have unique combination of pMHC and GEM.
    #    """
    #    self.plt_df = self.selected_mat.melt(ignore_index=False, var_name='gem', value_name='umi').dropna()
    #    self.plt_df.umi = self.plt_df.umi.astype(int)
        
    def transform_data_for_plotting(self):
        """
        For each row have unique combination of pMHC and GEM.
        """
        self.plt_df = self.selected_mat.melt(ignore_index=False, var_name='gem', value_name='umi').fillna(0) # OBS! new with replace!#.dropna()
        self.plt_df.umi = self.plt_df.umi.astype(int)

    def add_gem_count(self):
        """
        Make a variable for iterating over GEMs (fx in a GIF)
        """
        dct = (self.plt_df.groupby('gem', sort=False).size()
               .to_frame().reset_index().reset_index().set_index('gem')
               .rename(columns={'index':'gem_count',0:'query_count'}))
        dct['gem_count'] = dct.gem_count + 1

        # Make a var for iterating over GEMs
        self.plt_df['gem_count'] = self.plt_df.gem.map(dct.gem_count)

    def sort_data(self):
        self.plt_df.reset_index(inplace=True)
        if self.variable != 'sample_id':
            self.plt_df[self.var_lst] = self.plt_df[self.var_lst].astype("category") #why as category?
            self.plt_df[self.var_lst] = self.plt_df[self.var_lst].cat.set_categories(self.selected_queries)
        else:
            self.plt_df[self.var_lst] = self.plt_df[self.var_lst].fillna('').astype(str) #why do I have nans?
        self.plt_df.sort_values(by=self.var_lst, inplace=True)
            
    # Not needed anymore...
    #def add_scatter(self):
    #    for i, (name, subdf) in reversed(list(enumerate(self.plt_df.groupby('peptide_HLA_lst', sort=False)))): #reversed(list(grouped))
    #        self.plt_df.loc[subdf.index, 'y'] = np.random.normal(i, 0.1, subdf.shape[0])
    
    def transform_to_concordance(self):
        self.conc_df = (self.plt_df.sort_values(by=['gem','umi'])
                        .drop_duplicates(subset='gem', keep='last')
                        .groupby(self.var_lst).gem.count()
                        .to_frame())
        self.conc_df['clonotype'] = f'clonotype {self.ct}'
        self.conc_df.reset_index(inplace=True)
        self.conc_df['concordance'] = self.conc_df.gem / self.conc_df.gem.sum()
        self.conc_df.replace(0, np.nan, inplace=True)
         
    def plot_advanced_figure(self, figname):
        def get_legend_n_handle(l,h,key='gem'):
            leg = list()
            hdl = list()
            if key is None:
                keep = True
            else:
                keep = False
            for i,e in enumerate(l):
                if keep:
                    if int(float(e)) > 0:
                        leg.append(int(float(e)))
                        hdl.append(h[i])
                if e == key:
                    keep = True
            return hdl, leg
        
        fig = plt.figure(figsize=(20,7))
        fig.suptitle(f"Clonotype {self.ct}")

        gs = gridspec.GridSpec(1, 3, width_ratios=[2, 7, 2], wspace=0.2) #, left=0.05
        ax1 = fig.add_subplot(gs[0]) #ax1 = plt.subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1]) #ax2 = plt.subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[2]) #ax3 = plt.subplot(gs[0, 2])

        ########################
        # Add multipletplot
        ###############################
        tx1 = ax1.twinx()
        sns.scatterplot(data=self.plt_df, x="gem", y='peptide_HLA_lst', size='umi', color="gray", ax=ax1, legend='brief')
        
        h,l = ax1.get_legend_handles_labels()
        h,l = get_legend_n_handle(l, h, key=None)
        ax1.legend(h, l, bbox_to_anchor=(-1, 0.5), loc=10, frameon=False, title='UMI')
        
        sns.scatterplot(data=self.plt_df, x="gem", y='peptide_HLA_lst', size='umi', color="gray", ax=tx1, legend=False)

        ########################
        # Add boxplot
        ###############################
        PROPS = {'boxprops':{'alpha':0.3}, #'facecolor':'none', 
                 'medianprops':{'alpha':0.3},
                 'whiskerprops':{'alpha':0.3},
                 'capprops':{'alpha':0.3}}
        EMPTY = {'boxprops':{'alpha':0}, #'facecolor':'none', 
                 'medianprops':{'alpha':0},
                 'whiskerprops':{'alpha':0},
                 'capprops':{'alpha':0}}
        ARGS = {'x':"umi", 'y':"peptide_HLA_lst", 'data':self.plt_df, 'showfliers':False}
        order = self.plt_df.peptide_HLA_lst.unique()#[::-1]

        tx2 = ax2.twinx() # hack to get matching ticks on the right
        sns.boxplot(**ARGS, **PROPS, order=order, ax=ax2)
        sns.stripplot(data=self.plt_df, x="umi", y='peptide_HLA_lst', ax=ax2, order=order, jitter=0.2, edgecolor='white',linewidth=0.5, size=6)
        #sns.boxplot(**ARGS, **PROPS, ax=ax2, order=order) #, order=order
        sns.boxplot(**ARGS, **EMPTY, order=order, ax=tx2)
        #sns.scatterplot(data=self.plt_df, x="umi", y="y", hue='peptide_HLA_lst', legend=False, ax=ax2)

        # Add significance bar
        if self.test_dist():
            l = len(self.plt_df.peptide_HLA_lst.unique()) - 1
            y = [0,0,1,1] #[l, l, l-1, l-1] #
            x0 = self.plt_df.umi.max()
            x1 = x0 * 1.02
            x2 = x0 * 1.025
            x3 = x0 * 1.035

            ax2.plot([x1, x2, x2, x1], y, lw=0.7, c='0') #lw=1.5, 
            ax2.plot(x3, np.mean(y), marker="*", c='0')

        ######################################
        # Add concordance plot
        #########################################
        # Hack to get colorbar
        plot = ax3.scatter([np.nan]*len(order), order, c=[np.nan]*len(order), cmap='viridis_r', vmin=0, vmax=1) #[],[] #order[::-1]
        fig.colorbar(plot, ax=ax3)
        sns.scatterplot(data=self.conc_df, x='clonotype', y='peptide_HLA_lst',
                        size='gem', hue='concordance', hue_norm=(0,1), palette='viridis_r', ax=ax3) #g = 

        # Remove automatic sns legend for hue, keep only legend for size.
        h,l = ax3.get_legend_handles_labels()
        h,l = get_legend_n_handle(l, h, key='gem')
        ax3.legend(h, l, bbox_to_anchor=(1.5, 0.5), loc=6, frameon=False, title='GEM') #h[-4:], l[-4:]

        ######################################
        # Prettify
        #########################################

        ax1.set_title('Peptide MHC multiplets')
        ax2.set_title('UMI distribution per peptide MHC')
        ax3.set_title('Concordance')

        xmax = round(self.plt_df.umi.max(), -1)
        ax2.set_xticks(np.arange(0, xmax+10, 10))
        ax3.set_yticks([])

        ax1.set_xticklabels([])
        tx1.set_yticklabels([f'n= {n}' for n in self.plt_df.groupby('peptide_HLA_lst').size()])
        ax2.set_yticklabels([])
        ax3.set_yticklabels([])

        ax1.set_xlabel('GEM')
        ax1.set_ylabel('Peptide HLA', labelpad=20)
        tx1.set_ylabel('')
        ax2.set_xlabel('Peptide UMI')
        ax2.set_ylabel('')
        tx2.set_ylabel('')
        ax3.set_ylabel('')
        ax3.set_xlabel('')

        tx1.tick_params(axis='y', pad=20)

        ax2.spines['bottom'].set_bounds(0, xmax) # Hack to fix x-axis

        sns.despine(trim=True, right=True, ax=ax1)
        sns.despine(trim=True, right=False, ax=tx1)
        sns.despine(trim=True, right=False, ax=ax2)
        sns.despine(trim=True, right=False, ax=tx2)
        sns.despine(trim=True, left=True, ax=ax3) #offset={'left':-5,'right':-5}

        plt.savefig(figname %(self.fig_flg, self.ct), bbox_inches='tight')
        plt.show()
    
    
    def get_plotting_stats(self):
        self.summed_umis = self.summary_df.s.head(10)
        self.summed_gems = (self.mat > 0).sum(axis=1)
        
    def test_dist(self):
        assert not self.selected_mat.isna().any().any()
        # Select the peptides to test: the most abundant and the second most abundant (UMI wise)
        p1 = self.summary_df.index[0]
        p2 = self.summary_df.index[1]
        # Extract the UMI distribution for the two selected peptides.
        s1 = self.selected_mat.T[p1]#.fillna(0)
        s2 = self.selected_mat.T[p2]#.fillna(0)
        
        if sum(s1.fillna(0)-s2.fillna(0)) == 0:
            return False
        
        w, p = stats.wilcoxon(s1.fillna(0)-s2.fillna(0), alternative='greater')
        #u, p = stats.mannwhitneyu(s1, s2, alternative="greater") #, method="exact"
        
        if p <= 0.05:
            return True
        return False
        
    #def test_dist(self):
    #    """
    #    H0: No difference in UMI counts between the top two peptides
    #    H1: Top peptide has more UMI counts
    #
    #    If H1, then we assume that the clonotype has a clear preference for that peptide.
    #    """
    #    # Select the peptides to test: the most abundant and the second most abundant (UMI wise)
    #    p1 = self.summary_df.index[0]
    #    p2 = self.summary_df.index[1]
    #    # Extract the UMI distribution for the two selected peptides.
    #    s1 = self.mat.T[p1].dropna().to_list()
    #    s2 = self.mat.T[p2].dropna().to_list()
    #    # If fewer than 5 observations for any of the peptides, then assume H0
    #    n = min(len(s1), len(s2))
    #    if n < 5:
    #        return False
    #
    #    # Repeat experiment XX times to get more robust selection.
    #    p_list = list()
    #    for _ in range(50):
    #        l1 = sample(s1, n)
    #        l2 = sample(s2, n)
    #
    #        t, p = stats.ttest_ind(l1, l2, equal_var=False, alternative='greater')
    #
    #        p_list.append(p)
    #
    #    if np.median(p_list) <= 0.05:
    #        return True
    #    return False
    
    def get_imputed_query(self):
        return self.summary_df.index[0]
    
    def get_remaining_queries(self):
        return self.summary_df.index[1:10]
    
    def update_bins(self):
        Evaluate_Clonotype.value_bin.update([self.ct])
        Evaluate_Clonotype.trash_bin.update(self.summary_df.index[1:10])
        
    def update_variable_analysis(self):
        Evaluate_Clonotype.ct_checks[self.variable] = self.plt_df
        
    def update_flag(self, flag):
        self.fig_flg = flag

#########################################################################################################
#                                                 Input                                                 #
#########################################################################################################
INPUT = snakemake.input[0]
HASHING = snakemake.input[1]
OUTPUT = snakemake.output[0]
PLOT = snakemake.params[0] # Cannot use output, since I don't know which clonotypes will be significant
# PLOT = some_dir/%s/%d.pdf

########################################################################################################
#                                                 Load                                                 #
########################################################################################################
df = pd.read_csv(INPUT, converters=converters)
hsh = pd.read_excel(HASHING, sheet_name='HSH')

hsh.fillna('', inplace=True)
hsh['HLA'] = hsh['HLA_A'] + ', ' + hsh['HLA_B'] + ', ' + hsh['HLA_C']
hsh['HLA'] = hsh['HLA'].str.split(r',\s?').apply(lambda x: [i for i in x if i!= ''])
hsh.set_index('sample_id', inplace=True)

########################################################################################################
#                                               Prepare                                                #
########################################################################################################
if df.sample_id.nunique() == 0:
    variables = ['peptide_HLA']
    imp_vars = ['ct_pep']
    no_hashing = True
else:
    variables = ['sample_id','HLA_cd8','peptide_HLA']
    imp_vars = ['ct_sample','ct_hla','ct_pep']
    no_hashing = False

selected_clonotypes = df.groupby('ct').size()
selected_clonotypes = selected_clonotypes[selected_clonotypes >= 10] # OBS!
print('Number of clonotypes to test:', len(selected_clonotypes))
print(selected_clonotypes)

for ct, size in selected_clonotypes.items():
    print(ct, size)
    for variable, imp_var in zip(variables, imp_vars):
        #print(df[df.ct == ct].sample_id)
        #print(df[df.ct == ct].sample_id_lst)
        inst = Evaluate_Clonotype(df, ct, selected_clonotypes, variable=variable)
        inst.sum_umi()
        inst.calc_summary()
        #df.loc[instance.idx, 'umi_count_mhc_rel'] = instance.calc_relative_umi()
        inst.select_queries()
        inst.transform_data_for_plotting()
        inst.add_gem_count()
        inst.sort_data()
        inst.transform_to_concordance() # Only relevant for pMHC...
        inst.update_variable_analysis()

        if inst.test_dist():
            print(f'significant {variable}')
            df.loc[inst.idx, imp_var] = inst.get_imputed_query()
            
            if variable == 'peptide_HLA':
                inst.update_bins()
            
    if ct in Evaluate_Clonotype.value_bin:
        tmp = calc_binding_concordance(df[df.ct == ct].copy(), 'ct')
        conc_pep = tmp[tmp.binding_concordance == tmp.binding_concordance.max()].peptide_HLA.unique()[0]
        if conc_pep == tmp.ct_pep.unique()[0]:
            inst.update_flag('significant_match')
        else:
            inst.update_flag('significant_mismatch')
        print(f'{ct} is credible')
    else:
        inst.update_flag('insignificant')
       
    ########################################################################################################
    #                                                 Plot                                                 #
    ########################################################################################################
    #try:
    #if variable == 'peptide_HLA'
    #    inst.plot_advanced_figure('../lol_%s_%s.png') #../experiments/exp13/run2/plt/eval_clonotypes/significant/
    #except:
    #    print('not plottet')
    #finally:
    #    plt.close()
    #    plt.clf()
    
    if no_hashing:
        fig, ax1 = plt.subplots(1,1) #, figsize=(7,7)
        # pMHC
        order = Evaluate_Clonotype.ct_checks[variable].peptide_HLA_lst.unique() #inst.plt_df.peptide_HLA_lst.unique()#[::-1]
        ARGS = {'x':"umi", 'y':"peptide_HLA_lst", 'data':Evaluate_Clonotype.ct_checks['peptide_HLA'], 'order':order, 'ax':ax1}
        sns.boxplot(**ARGS, showfliers=False)
        sns.stripplot(**ARGS, jitter=0.2, edgecolor='white',linewidth=0.5, size=6)
        ax1.set_title("Peptide HLA")
        sns.despine(trim=True, ax=ax1)
        
        fig.suptitle(f'Clonotype {inst.ct} ({size} GEMs)')
        #plt.show()
        fig.savefig(PLOT %(inst.fig_flg, inst.ct), bbox_inches='tight')
        plt.close()
        plt.clf()
    else:
        ### PLOTTING ###
        fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(20,7))
        # pMHC
        order = Evaluate_Clonotype.ct_checks[variable].peptide_HLA_lst.unique() #inst.plt_df.peptide_HLA_lst.unique()#[::-1]
        ARGS = {'x':"umi", 'y':"peptide_HLA_lst", 'data':Evaluate_Clonotype.ct_checks['peptide_HLA'], 'order':order, 'ax':ax1}
        sns.boxplot(**ARGS, showfliers=False)
        sns.stripplot(**ARGS, jitter=0.2, edgecolor='white',linewidth=0.5, size=6)
        ax1.set_title("Peptide HLA")
        sns.despine(trim=True, ax=ax1)

        # Sample
        ARGS = {'x':"umi", 'y':"sample_id_lst", 'data':Evaluate_Clonotype.ct_checks['sample_id'], 'ax':ax3}
        sns.boxplot(**ARGS, showfliers=False)
        sns.stripplot(**ARGS, jitter=0.2, edgecolor='white',linewidth=0.5, size=6)
        ax3.set_title('Sample ID')
        sns.despine(trim=True, ax=ax3)
        #sns.scatterplot(**ARGS, x='gem')
        #sns.scatterplot(data=smpl.plt_df, x="gem", y='sample_id', size='umi', ax=ax1, legend='brief') #hue='sample_id', 

        # Hashing
        ARGS = {'x':"umi", 'y':"HLA_pool_cd8", 'data':Evaluate_Clonotype.ct_checks['HLA_cd8'], 'ax':ax2}
        sns.boxplot(**ARGS, showfliers=False)
        sns.stripplot(**ARGS, jitter=0.2, edgecolor='white',linewidth=0.5, size=6)
        ax2.set_title('Sample HLA')
        sns.despine(trim=True, ax=ax2)

        fig.suptitle(f'Clonotype {inst.ct} ({size} GEMs)')
        #plt.show()
        fig.savefig(PLOT %(inst.fig_flg, inst.ct), bbox_inches='tight')
        plt.close()
        plt.clf()

# Possible fix:
# https://github.com/mwaskom/seaborn/commit/1a537c100dd58c4a22187b8f2a02ab53a88030a2
# Check the sns version on computerome.
    
    
########################################################################################################
#                                                 Eval                                                 #
########################################################################################################
def notnan(var):
    return var == var

def determine_pep_match(row):
    if notnan(row.peptide_HLA) & notnan(row.ct_pep):
        return row.peptide_HLA == row.ct_pep
    else:
        return np.nan

def determine_hla_match_by_HLA_conc(row):
    if notnan(row.HLA_mhc) & (notnan(row.ct_hla) | notnan(row.sample_hla)):
        smp_hla = row.sample_hla if notnan(row.sample_hla) else []
        return (row.HLA_mhc == row.ct_hla) | (row.HLA_mhc in smp_hla) #row.sample_hla
    else:
        return np.nan

def determine_ct_match_by_HLA_conc(row):
    pep_hla = row.ct_pep.split(' ')[1] if notnan(row.ct_pep) else np.nan
    imp_hla = row.ct_hla
    smp_hla = row.sample_hla
    if notnan(pep_hla) & (notnan(imp_hla) | notnan(smp_hla)):
        smp_hla = smp_hla if notnan(smp_hla) else []
        return (pep_hla == imp_hla) | (pep_hla in smp_hla)
    else:
        return np.nan
    
def determine_train_label(row):
    if any((notnan(row.pep_match), notnan(row.hla_match), notnan(row.ct_match))):
        return (row.pep_match == True) & (row.hla_match == True) & (row.ct_match == True)
        #return (row.pep_match != False) & (row.hla_match != False) & (row.ct_match != False)
    else:
        return np.nan
    
def determine_test_label(row):
    if notnan(row.HLA_mhc) & notnan(row.HLA_cd8):
        return row.HLA_mhc in row.HLA_cd8
    else:
        return np.nan


try:
    df['sample_hla'] = df.ct_sample.map(hsh.HLA)

    df['pep_match'] = df.apply(lambda row: determine_pep_match(row), axis=1) #df.ct_pep == df.peptide_HLA
    df['hla_match'] = df.apply(lambda row: determine_hla_match_by_HLA_conc(row), axis=1)
    df['ct_match'] = df.apply(lambda row: determine_ct_match_by_HLA_conc(row), axis=1)

    # A clonotype is only credible when we significantly can determine a peptide and an HLA
    # This might be a harsh criteria for some datasets?!
    df['valid_ct'] = df.ct.isin(Evaluate_Clonotype.value_bin) & (df.ct_match == True) #(df.ct_match != False)

    # A TP within the trainset requires a significance and concordance of pMHC and sHLA
    df['train_label'] = df.apply(lambda row: determine_train_label(row), axis=1)
    df['test_label'] = df.apply(lambda row: determine_test_label(row), axis=1) 
except AttributeError:
    df['sample_hla'] = np.nan

    df['pep_match'] = df.apply(lambda row: determine_pep_match(row), axis=1) #df.ct_pep == df.peptide_HLA
    df['hla_match'] = np.nan
    df['ct_match'] = np.nan

    # A clonotype is only credible when we significantly can determine a peptide and an HLA
    # This might be a harsh criteria for some datasets?!
    df['valid_ct'] = df.ct.isin(Evaluate_Clonotype.value_bin) #& (df.ct_match == True) #(df.ct_match != False)

    # A TP within the trainset requires a significance and concordance of pMHC and sHLA
    df['train_label'] = df.apply(lambda row: determine_train_label(row), axis=1)
    df['test_label'] = df.apply(lambda row: determine_test_label(row), axis=1)
df.to_csv(OUTPUT, index=False)

#valid_df = df[df.ct.isin(Evaluate_Clonotype.value_bin)].fillna(0)
#valid_df.to_csv(OUTPUT, index=False)
#df.to_csv('/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run2/cat/tables/tcr_barcode.cleaned.valid_ct.csv', index=False)

