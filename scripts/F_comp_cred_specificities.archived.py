#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from ast import literal_eval
import re
from scipy import stats
from random import sample

plt.style.use('ggplot')

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
        self.ct = ct
        self.idx = grp.index
        self.rel_umi = use_relative_umi
        self.fig_flg = 'None'
        self.fig, (self.ax1, self.ax2) = plt.subplots(1,2,figsize=(10,3))
        
        # Initialize list of clonotypes
        self.sel_cts = selected_clonotypes
        self.sel_cts.index = self.sel_cts.index.astype(int)
        
        # Initialize matrix of pMHC count per GEM (not method (df is not self))
        self.all_peptide_hla = df.peptide_HLA_lst.explode().unique()
        self.mat = pd.DataFrame(index=all_peptide_HLA, columns=df.gem.unique())
        
        # Count no. of GEMs within grp that are annotated with a specific pMHC (method?)
        self.gems_per_pep_ann = grp.groupby('peptide_HLA').size() # annotated_gems
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
            self.mat = self.mat / self.grp.umi_count_mhc.quantile(0.9, interpolation='lower')
        return self.grp.umi_count_mhc / self.grp.umi_count_mhc.quantile(0.9, interpolation='lower')
    
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

        ## Alternative to subsampling (not implemented)
        #_, pp = stats.ttest_ind(s1, s2, equal_var=False, alternative='greater')
        #
        #if pp <= 0.05:
        #    print('*** Significant with no sampling')
        #else:
        #    print('--- Insignificant even with no sampling')
    
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
                           %(self.ct, self.selected_clonotypes[self.ct],
                             self.grp.genes_TRA.unique()[0], self.grp.cdr3_TRA.unique()[0],
                             self.grp.genes_TRB.unique()[0], self.grp.cdr3_TRB.unique()[0]))

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

for ct in selected_clonotypes:
    instance = Evaluate_Clonotype(df, ct, selected_clonotypes, use_relative_umi=False)
    instance.sum_umi()
    instance.calc_summary()
    
    df.loc[instance.idx, 'umi_count_mhc_rel'] = instance.calc_relative_umi()
    
    instance.get_plotting_stats()
    
    # test if two annotations are significantly different distributed.
    if instance.test_dist():
        df.loc[instance.idx, 'ct_pep'] = instance.get_imputed_peptide()
        
        instance.update_bins()
        instance.update_flag('significant')
    else:
        instance.update_flag('insignificant')
        
    instance.plot_histogram()
    instance.plot_boxplot()
    instance.add_labels()
    instance.save_figure(PLOT)

df['pep_match'] = df.ct_pep == df.peptide_HLA

valid_df = df[df.ct.isin(Evaluate_Clonotype.value_bin)].fillna(0)
valid_df.to_csv(OUTPUT, index=False)
