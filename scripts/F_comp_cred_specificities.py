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
import argparse

import os
import sys

sns.set_style('ticks', {'axes.edgecolor': '0',  
                        'xtick.color': '0',
                        'ytick.color': '0'})
sns.set_context("paper",font_scale=2)

def lst_converter(x):
    return x.split("|")

def HLA_lst_converter(x):
    return [hla_lst.split(';') for hla_lst in x.split('|')]

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
    try:
        return [] if val == '' else [v for v in literal_eval(val)]
    except:
        return [] if val == '' else [float(v) for v in lst_converter(val)]
    

converters = {'peptide_HLA_lst': lst_converter, #peptide_hla_converter,
              'umi_count_lst_mhc': literal_converter, #literal_eval,
              'umi_count_lst_cd8': literal_converter,
              'umi_count_lst_TRA': literal_converter,'umi_count_lst_TRB': literal_converter,
              'cdr3_lst_TRA': lst_converter, #cdr3_lst_converter,
              'cdr3_lst_TRB': lst_converter, #cdr3_lst_converter,
              'HLA_lst_mhc': lst_converter, #cdr3_lst_converter,
              'HLA_pool_cd8': lst_converter, #cdr3_lst_converter,
              'HLA_cd8': lst_converter, #HLA_cd8_converter,
              'HLA_lst_cd8': HLA_lst_converter,'sample_id_lst':literal_converter}

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

def get_argparser():
    """Return an argparse argument parser."""
    parser = argparse.ArgumentParser(prog = 'Evaluate clonotypes',
                                     description = 'Evaluates which clonotypes can be used for grid search to identify optimal UMI thresholds')
    add_arguments(parser)
    return parser

def add_arguments(parser):
    parser.add_argument('--input', required=True, help='Filepath for data')
    parser.add_argument('--output', required=True, help='Filepath for output data')
    parser.add_argument('--plots', required=True, help='Filepath for output plots (must contain two placeholders, ie plt_dir/%s/%d.pdf)')

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
        
        # Initialize matrix of query count per GEM 
        self.queries = df[self.var_lst].explode().drop_duplicates()
        self.mat = pd.DataFrame(index=self.queries, columns=df.gem.unique())
        
        # Count no. of GEMs within grp that are annotated with a specific pMHC
        self.gems_per_query = self.df.explode(self.variable).groupby(self.variable).size()
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
        self.mat.fillna(0, inplace=True)
                
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
        
    def select_queries(self, n=11):
        self.selected_queries = self.summary_df.head(n).index
        self.selected_mat = self.mat.loc[self.mat.index.isin(self.selected_queries), self.df.gem]
        
    def transform_data_for_plotting(self):
        """
        For each row have unique combination of pMHC and GEM.
        """
        self.plt_df = self.selected_mat.melt(ignore_index=False, var_name='gem', value_name='umi').fillna(0)
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
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax3 = fig.add_subplot(gs[2])

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
        PROPS = {'boxprops':{'alpha':0.3},
                 'medianprops':{'alpha':0.3},
                 'whiskerprops':{'alpha':0.3},
                 'capprops':{'alpha':0.3}}
        EMPTY = {'boxprops':{'alpha':0},
                 'medianprops':{'alpha':0},
                 'whiskerprops':{'alpha':0},
                 'capprops':{'alpha':0}}
        ARGS = {'x':"umi", 'y':"peptide_HLA_lst", 'data':self.plt_df, 'showfliers':False}
        order = self.plt_df.peptide_HLA_lst.unique()#[::-1]

        tx2 = ax2.twinx() # hack to get matching ticks on the right
        sns.boxplot(**ARGS, **PROPS, order=order, ax=ax2)
        sns.stripplot(data=self.plt_df, x="umi", y='peptide_HLA_lst', ax=ax2, order=order, jitter=0.2, edgecolor='white',linewidth=0.5, size=6)
        sns.boxplot(**ARGS, **EMPTY, order=order, ax=tx2)

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
        plot = ax3.scatter([np.nan]*len(order), order, c=[np.nan]*len(order), cmap='viridis_r', vmin=0, vmax=1)
        fig.colorbar(plot, ax=ax3)
        sns.scatterplot(data=self.conc_df, x='clonotype', y='peptide_HLA_lst',
                        size='gem', hue='concordance', hue_norm=(0,1), palette='viridis_r', ax=ax3)

        # Remove automatic sns legend for hue, keep only legend for size.
        h,l = ax3.get_legend_handles_labels()
        h,l = get_legend_n_handle(l, h, key='gem')
        ax3.legend(h, l, bbox_to_anchor=(1.5, 0.5), loc=6, frameon=False, title='GEM')

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
        sns.despine(trim=True, left=True, ax=ax3)

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
        s1 = self.selected_mat.T[p1]
        s2 = self.selected_mat.T[p2]
        
        if sum(s1.fillna(0)-s2.fillna(0)) == 0:
            return False
        
        w, p = stats.wilcoxon(s1.fillna(0)-s2.fillna(0), alternative='greater')
        
        if p <= 0.05:
            return True
        return False
    
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
try:
    INPUT = snakemake.input.data
    OUTPUT = snakemake.output.data
    PLOT = snakemake.params.plots
except:
    parser = get_argparser()
    args = parser.parse_args()

    INPUT = args.input
    OUTPUT = args.output
    PLOT = args.plots

########################################################################################################
#                                                 Load                                                 #
########################################################################################################
df = pd.read_csv(INPUT, converters=converters)

########################################################################################################
#                                               Prepare                                                #
########################################################################################################

variables = ['peptide_HLA']
imp_vars = ['ct_pep']
no_hashing = True

selected_clonotypes = df.groupby('ct').size()
selected_clonotypes = selected_clonotypes[selected_clonotypes >= 10] # OBS!

for ct, size in selected_clonotypes.items():
    for variable, imp_var in zip(variables, imp_vars):
        inst = Evaluate_Clonotype(df, ct, selected_clonotypes, variable=variable)
        inst.sum_umi()
        inst.calc_summary()
        inst.select_queries()
        inst.transform_data_for_plotting()
        inst.add_gem_count()
        inst.sort_data()
        inst.transform_to_concordance() # Only relevant for pMHC...
        inst.update_variable_analysis()

        if inst.test_dist():
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
    else:
        inst.update_flag('insignificant')
       
    ########################################################################################################
    #                                                 Plot                                                 #
    ########################################################################################################
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

        # Hashing
        ARGS = {'x':"umi", 'y':"HLA_pool_cd8", 'data':Evaluate_Clonotype.ct_checks['HLA_cd8'], 'ax':ax2}
        sns.boxplot(**ARGS, showfliers=False)
        sns.stripplot(**ARGS, jitter=0.2, edgecolor='white',linewidth=0.5, size=6)
        ax2.set_title('Sample HLA')
        sns.despine(trim=True, ax=ax2)

        fig.suptitle(f'Clonotype {inst.ct} ({size} GEMs)')
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


df['HLA_mhc'] = df.peptide_HLA.str.split(' ', expand=True)[1]
df['pep_match'] = df.apply(lambda row: determine_pep_match(row), axis=1)
df['valid_ct'] = df.ct.isin(Evaluate_Clonotype.value_bin)
    
df.to_csv(OUTPUT, index=False)

