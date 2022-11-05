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

import os
import sys
from D_plot_specificity_matrix_utils import (calc_binding_concordance)

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
              'HLA_lst_mhc': cdr3_lst_converter,'HLA_cd8': HLA_cd8_converter}

#########################################################################################################
#                                                 Class                                                 #
#########################################################################################################
class Evaluate_Clonotype():
    """
    An instance is a clonotype subset of the data.
    """
    value_bin = set() # clonotypes
    trash_bin = set() # peptide-HLAs
    
    def __init__(self, df, ct, selected_clonotypes, use_relative_umi=False):
        self.df = df[df.ct == ct]
        self.ct = int(ct)
        self.idx = self.df.index
        self.rel_umi = use_relative_umi
        self.fig_flg = 'None'
        #self.fig, (self.ax1, self.ax2) = plt.subplots(1,2,figsize=(10,3))
        
        # Initialize list of clonotypes
        self.sel_cts = selected_clonotypes
        self.sel_cts.index = self.sel_cts.index.astype(int)
        
        # Initialize matrix of pMHC count per GEM (not method (df is not self))
        self.all_peptide_hla = df.peptide_HLA_lst.explode().drop_duplicates() #.unique()
        self.mat = pd.DataFrame(index=self.all_peptide_hla, columns=df.gem.unique())
        
        # Count no. of GEMs within grp that are annotated with a specific pMHC (method?)
        self.gems_per_pep_ann = self.df.groupby('peptide_HLA').size() # annotated_gems
        self.gems_per_pep_all = pd.concat([self.gems_per_pep_ann,
                                           pd.Series(0, index=self.all_peptide_hla[~self.all_peptide_hla.isin(self.gems_per_pep_ann.index)])])
        
    def sum_umi(self):
        # Sum UMIs for each peptide across GEMs (multiple peptides per GEM)
        for idx, row in self.df.iterrows():
            if len(row.umi_count_lst_mhc) == len(row.peptide_HLA_lst):
                self.mat.loc[row.peptide_HLA_lst, row.gem] = row.umi_count_lst_mhc
            else:
                self.mat.loc[row.peptide_HLA_lst, row.gem] = [0] * len(row.peptide_HLA_lst)
                
    def calc_summary(self):
        self.summary_df = self.mat.sum(axis=1).sort_values(ascending=False).to_frame().rename(columns={0:'s'})
        self.summary_df['avg'] = self.mat.mean(axis=1)
        self.summary_df['col'] = 'grey'
        self.summary_df['r'] = self.summary_df.s / self.summary_df.s.max() # Unnecessary
        
    def calc_relative_umi(self):
        if self.rel_umi:
            self.mat = self.mat / self.df.umi_count_mhc.quantile(0.9, interpolation='lower')
        return self.df.umi_count_mhc / self.df.umi_count_mhc.quantile(0.9, interpolation='lower')
    
    def select_peptides(self, n=11):
        self.selected_peptides = self.summary_df.head(n).index
        self.selected_mat = self.mat[self.mat.index.isin(self.selected_peptides)]

    def transform_data_for_plotting(self):
        """
        For each row have unique combination of pMHC and GEM.
        """
        self.plt_df = self.selected_mat.melt(ignore_index=False, var_name='gem', value_name='umi').dropna()

    def add_gem_count(self):
        """
        Make a variable for iterating over GEMs (fx in a GIF)
        """
        dct = (self.plt_df.groupby('gem', sort=False).size()
               .to_frame().reset_index().reset_index().set_index('gem')
               .rename(columns={'index':'gem_count',0:'peptide_count'}))
        dct['gem_count'] = dct.gem_count + 1

        # Make a var for iterating over GEMs
        self.plt_df['gem_count'] = self.plt_df.gem.map(dct.gem_count)

    def sort_data(self):
        self.plt_df.reset_index(inplace=True)
        self.plt_df.peptide_HLA_lst = self.plt_df.peptide_HLA_lst.astype("category")
        self.plt_df['peptide_HLA_lst'] = self.plt_df.peptide_HLA_lst.cat.set_categories(self.selected_peptides) # instead of inplace
        self.plt_df.sort_values(by='peptide_HLA_lst', inplace=True)

    def add_scatter(self):
        for i, (name, subdf) in enumerate(self.plt_df.groupby('peptide_HLA_lst', sort=False)): #reversed(list(grouped))
            self.plt_df.loc[subdf.index, 'y'] = np.random.normal(i, 0.1, subdf.shape[0])
    
    def transform_to_concordance(self):
        self.conc_df = (self.plt_df.sort_values(by=['gem','umi'])
                        .drop_duplicates(subset='gem', keep='last')
                        .groupby('peptide_HLA_lst').gem.count()
                        .to_frame())
        self.conc_df['clonotype'] = f'clonotype {self.ct}'
        self.conc_df.reset_index(inplace=True)
        self.conc_df['concordance'] = self.conc_df.gem / self.conc_df.gem.sum()
        self.conc_df.replace(0, np.nan, inplace=True)
         
    def plot_advanced_figure(self, figname):
        fig = plt.figure(figsize=(20,7))
        fig.suptitle(f"Clonotype {ct}")

        gs = gridspec.GridSpec(1, 3, width_ratios=[2, 7, 2], wspace=0.2) #, left=0.05
        ax1 = fig.add_subplot(gs[0]) #ax1 = plt.subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1]) #ax2 = plt.subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[2]) #ax3 = plt.subplot(gs[0, 2])

        ########################
        # Add multipletplot
        ###############################
        tx1 = ax1.twinx()
        sns.scatterplot(data=self.plt_df, x="gem", y='peptide_HLA_lst', size='umi', color="gray", ax=ax1)
        sns.scatterplot(data=self.plt_df, x="gem", y='peptide_HLA_lst', size='umi', color="gray", ax=tx1, legend=False)
        ax1.legend(bbox_to_anchor=(-1, 0.5), loc=10, frameon=False, title='UMI')

        ###############################
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

        tx2 = ax2.twinx() # hack to get matching ticks on the right
        sns.boxplot(**ARGS, **PROPS, ax=ax2)
        sns.boxplot(**ARGS, **EMPTY, ax=tx2)
        sns.scatterplot(data=self.plt_df, x="umi", y="y", hue='peptide_HLA_lst', legend=False, ax=ax2)

        # Add significance bar
        if self.test_dist():
            y = [0,0,1,1]
            x0 = self.plt_df.umi.max()
            x1 = x0 * 1.02
            x2 = x0 * 1.025
            x3 = x0 * 1.035

            ax2.plot([x1, x2, x2, x1], y, lw=0.7, c='0') #lw=1.5, 
            ax2.plot(x3, np.mean(y), marker="*", c='0')

        #########################################
        # Add concordance plot
        #########################################
        # Hack to get colorbar
        plot = ax3.scatter([], [], c=[], cmap='viridis_r', vmin=0, vmax=1)
        fig.colorbar(plot, ax=ax3)
        sns.scatterplot(data=self.conc_df, x='clonotype', y='peptide_HLA_lst',
                        size='gem', hue='concordance', hue_norm=(0,1), palette='viridis_r', ax=ax3) #g = 

        # Remove automatic sns legend for hue, keep only legend for size.
        h,l = ax3.get_legend_handles_labels()
        ax3.legend(h[-4:], l[-4:], bbox_to_anchor=(1.5, 0.5), loc=6, frameon=False, title='GEM')

        #########################################
        # Prettify
        #########################################

        ax1.set_title('Peptide MHC multiplets')
        ax2.set_title('UMI distribution per peptide MHC')
        ax3.set_title('Concordance')

        xmax = round(self.plt_df.umi.max(), -1)
        if xmax > 250:
            interval = 20
        else:
            interval = 10
        ax2.set_xticks(np.arange(0, xmax+interval, interval))
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

        plt.savefig(name %(self.fig_flg, self.ct), bbox_inches='tight')
        #plt.show()
    
    
    def get_plotting_stats(self):
        self.summed_umis = self.summary_df.s.head(10)
        self.summed_gems = (self.mat > 0).sum(axis=1)
        
    def test_dist(self):
        """
        H0: No difference in UMI counts between the top two peptides
        H1: Top peptide has more UMI counts

        If H1, then we assume that the clonotype has a clear preference for that peptide.
        """
        # Select the peptides to test: the most abundant and the second most abundant (UMI wise)
        p1 = self.summary_df.index[0]
        p2 = self.summary_df.index[1]
        # Extract the UMI distribution for the two selected peptides.
        s1 = self.mat.T[p1].dropna().to_list()
        s2 = self.mat.T[p2].dropna().to_list()
        # If fewer than 5 observations for any of the peptides, then assume H0
        n = min(len(s1), len(s2))
        if n < 5:
            return False

        # Repeat experiment XX times to get more robust selection.
        p_list = list()
        for _ in range(50):
            l1 = sample(s1, n)
            l2 = sample(s2, n)

            t, p = stats.ttest_ind(l1, l2, equal_var=False, alternative='greater')

            p_list.append(p)

        if np.median(p_list) <= 0.05:
            return True
        return False
    
    def get_imputed_peptide(self):
        return self.summary_df.index[0]
    
    def get_remaining_peptides(self):
        return self.summary_df.index[1:10]
    
    def update_bins(self):
        Evaluate_Clonotype.value_bin.update([self.ct])
        Evaluate_Clonotype.trash_bin.update(self.summary_df.index[1:10])
        
    def update_flag(self, flag):
        self.fig_flg = flag
        
    def plot_histogram(self):
        for peptide_HLA in self.summed_umis.index:
            data = self.mat.T[peptide_HLA].dropna()
            labl = '%s %s %d (%d) (%d)' %(peptide_HLA.split(' ')[1],
                                          peptide_HLA.split(' ')[0],
                                          self.summed_umis[peptide_HLA],
                                          self.summed_gems[peptide_HLA],
                                          self.gems_per_pep_all[peptide_HLA])
            # Plot
            self.ax2.hist(data, alpha=0.2, label=labl)
        # Save plotting colors
        color_idx = self.summary_df.columns.to_list().index('col')
        print(len(plt.rcParams['axes.prop_cycle'].by_key()['color']))
        print(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        self.summary_df.iloc[0:10, color_idx] = plt.rcParams['axes.prop_cycle'].by_key()['color']
        
    def plot_boxplot(self):
        boxes = list()
        self.meaned_umis = self.summary_df.sort_values(by='avg',ascending=False).head(10)
        for peptide_HLA in self.meaned_umis.index:
            boxes.append(self.mat.T[peptide_HLA].dropna().to_list())
        bp = self.ax1.boxplot(boxes, patch_artist=True) #, vert=False
        colors = self.meaned_umis.col #.to_list()

        # whiskers and caps have to be treated separately since there are two of each for each plot
        for item in ['whiskers', 'caps']:
            for sub_item, color in zip(zip(bp[item][::2], bp[item][1::2]), colors):
                plt.setp(sub_item, color=color, alpha=0.2)
        for item in ['fliers']:
            for sub_item, color in zip(bp[item], colors):
                plt.setp(sub_item, markeredgecolor=color, alpha=0.2)
        for item in ['boxes', 'medians']:
            for sub_item, color in zip(bp[item], colors):
                plt.setp(sub_item, color=color, alpha=0.2)

    def add_labels(self):
        # Boxplot
        self.ax1.set_title('Sorted by mean of BC UMI')
        self.ax1.set_xticklabels(self.meaned_umis.index, rotation=90)
        self.ax1.set_ylabel('UMI count')
        self.ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)
        self.ax1.set_axisbelow(True) # Hide these grid behind plot objects
        
        plt.legend(bbox_to_anchor=(1.05, 0.5),
                   loc='center left',
                   title='HLA peptide UMI-sum (found in #GEMs) (ann. in #GEMs)')
        
        # Histogram
        self.ax2.set_xlabel('UMI count')
        self.ax2.set_ylabel('GEM count')
        self.ax2.set_title('Clonotype %d (%d GEMs)\n%s - %s\n%s - %s'
                           %(self.ct, self.sel_cts[self.ct],
                             self.df.genes_TRA.unique()[0], self.df.cdr3_TRA.unique()[0],
                             self.df.genes_TRB.unique()[0], self.df.cdr3_TRB.unique()[0]))

    def save_figure(self, name):
        plt.savefig(name %(self.fig_flg, self.ct), bbox_inches='tight')

#########################################################################################################
#                                                 Input                                                 #
#########################################################################################################
INPUT = snakemake.input[0]
OUTPUT = snakemake.output[0]
PLOT = snakemake.params[0] # Cannot use output, since I don't know which clonotypes will be significant
# PLOT = some_dir/%s/%d.pdf

########################################################################################################
#                                                 Load                                                 #
########################################################################################################
df = pd.read_csv(INPUT, converters=converters)

########################################################################################################
#                                                 Main                                                 #
########################################################################################################
selected_clonotypes = df.groupby('ct').size()

for ct in selected_clonotypes.index:
    instance = Evaluate_Clonotype(df, ct, selected_clonotypes, use_relative_umi=False)
    instance.sum_umi()
    instance.calc_summary()
    
    df.loc[instance.idx, 'umi_count_mhc_rel'] = instance.calc_relative_umi()

    if instance.test_dist():
        df.loc[instance.idx, 'ct_pep'] = instance.get_imputed_peptide()

        instance.update_bins()
        instance.update_flag('significant')
    else:
        instance.update_flag('insignificant')
    
    instance.select_peptides()
    instance.transform_data_for_plotting()
    instance.add_gem_count()
    instance.sort_data()
    instance.add_scatter()
    instance.transform_to_concordance()
    instance.plot_advanced_figure(PLOT)

df['pep_match'] = df.ct_pep == df.peptide_HLA

valid_df = df[df.ct.isin(Evaluate_Clonotype.value_bin)].fillna(0)
valid_df.to_csv(OUTPUT, index=False)

