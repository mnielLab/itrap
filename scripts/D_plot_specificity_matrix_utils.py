#!/usr/bin/env python
# coding: utf-8

# In[1]:

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import itertools
import collections

plt.style.use('ggplot')

# # Functions
def epitope_sorter_index(df):
    EPITOPE_SORTER = ['CLYBL',
                      'v9', 'v15', 'v19', 'v3', 'v5', 'v6', 'v10', 'v13',
                      'v16', 'v24', 'v41', 'v2', 'v11', 'v18', 'v23', 'v25',
                      'v26', 'v27', 'v4', 'v7', 'v8', 'v12', 'v14', 'v20',
                      'v21', 'v1', 'v17', 'v22' ,'v30', 'v31', 'v36', 'v37',
                      'v38', 'v39', 'v40', 'v32', 'v33', 'v34', 'v35']

    sorterIndex = dict(zip(EPITOPE_SORTER,range(len(EPITOPE_SORTER))))
    
    return df.epitope.map(sorterIndex) #df['epitope_rank'] = 

def filtering(df, read_threshold, tcr_threshold,
              exclude_clonotype_singlets=False,
              exclude_specificity_singlets=False,
              clonotype_fmt='ct',
              filtration='no_filtration',
              umi_delta=0):

    df = df[df[filtration] & (df.delta_umi_mhc >= umi_delta)] #OBS! should alpha and beta also be filtered on UMI? & (df.delta_umi_TRA > umi_delta)

    if df.empty:
        print('empty df')
        return df

    if (filtration == 'no_filtration'):
        def f(row):
            if row.umi_count_mhc >= read_threshold:
                if (row.umi_count_TRA >= tcr_threshold) | (row.umi_count_TRB >= tcr_threshold):
                    return True
            return False
        flt_df = df[df.apply(lambda row: f(row), axis=1).values]
    else:
        flt_df = df[(df.umi_count_mhc >= read_threshold) &
                    (df.umi_count_TRA >= tcr_threshold) &
                    (df.umi_count_TRB >= tcr_threshold)]

    if exclude_clonotype_singlets:
        flt_df = flt_df[flt_df.duplicated(subset=clonotype_fmt, keep=False)]
    if exclude_specificity_singlets:
        flt_df = flt_df[flt_df.duplicated(subset=[clonotype_fmt, 'peptide_HLA'], keep=False)]
    return flt_df

def calc_binding_concordance(df, clonotype_fmt):
    #assert df.size > 0, "df empty"
    gems_per_specificity        = df.groupby([clonotype_fmt,'peptide_HLA']).gem.count().to_dict()
    df['gems_per_specificity']  = df.set_index([clonotype_fmt,'peptide_HLA']).index.map(gems_per_specificity)
    gems_per_spec_hla_match     = df[df.HLA_match == True].groupby([clonotype_fmt, 'peptide_HLA']).gem.count().to_dict()
    df['gems_per_spec_hla_match'] = df.set_index([clonotype_fmt,'peptide_HLA']).index.map(gems_per_spec_hla_match)
    gems_per_clonotype          = df.groupby([clonotype_fmt]).gem.count().to_dict()
    df['gems_per_clonotype']    = df[clonotype_fmt].map(gems_per_clonotype)
    df['binding_concordance']   = df.gems_per_specificity / df.gems_per_clonotype
    df['hla_concordance']       = df.gems_per_spec_hla_match / df.gems_per_specificity
    #print('any 0 gems_per_spec_hla_match?', any(df.gems_per_spec_hla_match == 0))
    #print('any NaN gems_per_spec_hla_match?', any(df.gems_per_spec_hla_match.isna()))
    #print('any NaN hla_concordance?', any(df.hla_concordance.isna()))
    df['hla_concordance']       = df.hla_concordance.fillna(0)
    return df

def get_nonsinglet_idxs(duplicated_list):
    dups = collections.defaultdict(list)
    for i, e in enumerate(duplicated_list):
        dups[e].append(i)
    idxs = []
    for k, v in dups.items():
        if len(v) > 1:
            idxs += v
    return idxs

class Plotting:
    def __init__(self, bc_threshold, tcr_threshold, plot_element, gems, tcrs, point_scale=1, grid_axis='both'):
        import matplotlib as mpl
        mpl.rcParams['axes.grid.axis'] = grid_axis #'both'#y
        self.fig, self.ax = plt.subplots(figsize=(20, 10))
        self.bc_threshold = bc_threshold
        self.tcr_threshold = tcr_threshold
        self.plot_element = plot_element
        self.unique_tcrs = tcrs.unique()
        self.unique_gems = gems.unique()
        self.point_scale = point_scale

    def var2col(self, var):
        col = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2']
        self.var2col = dict(zip(var.unique(), col*len(var.unique())))

    def get_color(self, var):
        return list(map(self.var2col.get, var))

    def add_labels_and_ticks(self, by='TCR', ylabel="Peptide, HLA"):
        #self.unique_tcrs = unique_tcrs
        #self.unique_gems = unique_gems

        if by == 'TCR':
            plt.tick_params(labelsize=8) 
            plt.xticks(rotation=90, size=2)
            plt.xlabel("%i Clonotypes (of %i GEMs)" %(len(self.unique_tcrs), len(self.unique_gems)), fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            plt.title("Specificity per clonotype", fontsize=20)
            plt.xlim(-1, len(self.unique_tcrs))
            plt.ylim(-1, len(self.ax.get_yticks()))
        elif by == 'GEM':
            plt.tick_params(labelbottom=False, labelright=True, labelsize=8)
            plt.xlabel("%i GEMs (of %i clonotypes)" %(len(self.unique_gems), len(self.unique_tcrs)), fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            plt.title("Specificity per GEM", fontsize=20)
            plt.xlim(-1, len(self.unique_gems))
            plt.ylim(-1, len(self.ax.get_yticks()))
        else:
            plt.tick_params(labelbottom=True, labelright=True, labelsize=8)
            plt.xticks(rotation=90)
            plt.xlabel("%i GEMs (of %i clonotypes)" %(len(self.unique_gems), len(self.unique_tcrs)), fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            plt.title("UMI count distributions", fontsize=20)
            #if self.df.umi_count_tcr.max() == self.df.umi_count_tcr.max(): # Test if max is NaN
            #    plt.ylim(0, int(self.df.umi_count_tcr.max())+1)

    def add_colorbar(self):
        sm = plt.cm.ScalarMappable(cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1))
        if mpl.__version__ < '3.1':
            sm._A = []
        plt.colorbar(sm)
        #plt.ylabel('Binding concordance', rotation=270)

    def add_filtering_legend(self):
        # Criteria
        textstr = '\n'.join((
            "Criteria",
            "Min. BC UMI count\t%i".expandtabs() %self.bc_threshold,
            "Min. TCR UMI count\t%i".expandtabs() %self.tcr_threshold))#,
            #"Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
            #"Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
        props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
        self.ax.text(0.05, 0.95,textstr,
            fontsize=12,
            horizontalalignment='left',
            verticalalignment='top',
            clip_on=False,
            transform=self.ax.transAxes,
            bbox=props)

    def add_label_legend(self, columns=15, xanchor=0.01):
        from matplotlib.lines import Line2D

        legend_elements = []
        for var, color in self.var2col.items():
            legend_elements += [Line2D([0],[0],
                                marker='o', color='w',
                                label=var,
                                markerfacecolor=color,
                                markersize=10)]

        legend1 = self.ax.legend(handles=legend_elements, ncol=columns, loc=2, bbox_to_anchor=(xanchor, -0.03)) #0.02
        self.ax.add_artist(legend1)

    def add_size_legend(self):
        if self.plot_element == 'GEMs':
            sizes = [1, 5, 10]
        else:
            sizes = [1, 2, 3]
        for size in sizes:
            plt.scatter([], [], c='k', alpha=0.3, s=size*self.point_scale, label='{} {}'.format(size, self.plot_element))
        plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title=self.plot_element, loc='lower right') #'upper left'

        #for marker, response in [('+', "TCR singlet")]: #('o', True), 
        #    plt.scatter([], [], c='k', alpha=0.3, s=50, label=str(response), marker=marker)

    def finalize(self, show, save_tuba, save_sund):
        if show:
            #plt.savefig(show + 'specificity_matrix.pdf', bbox_inches='tight')
            #np.savetxt(show + 'unique_gems.lst', list(SpecificityMatrix.unique_gems), fmt='%s')
            plt.show()
            print("OBS! Figures are not saved!")
            return
        if save_tuba:
            plt.savefig(save_tuba, bbox_inches='tight')
        if save_sund:
            plt.savefig(save_sund + project + clonotype_fmt + "/%s/"%filtration + "b%i.t%i.ecs_%s.ess_%s.pdf" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
        plt.cla()   # Clear axis
        self.fig.clf()   # Clear figure
        plt.close(self.fig) # Close a figure window
        #plt.cla()   # Clear axis
        #plt.clf()   # Clear figure
        #plt.close(self.fig) # Close a figure window




class SpecificityMatrix:
    remaining_indexes = set()

    def __init__(self, df, sub_df, mhc_barcode, clonotype_fmt):
        self.df = df
        self.sub_df = sub_df
        self.mhc_barcode = mhc_barcode
        self.clonotype_fmt = clonotype_fmt
        self.sub_clonotype = sub_df[clonotype_fmt]

    @property
    def relevant_indexes(self):
        """
        This method ensures that a group of GEMs is plotted within
        the plateau if it has the highest concordance of all
        the groups of GEMs belonging to the same clonotype
        """
        #print(self.sub_df)
        all_concordance_scores = self.df[self.df[self.clonotype_fmt].isin(self.sub_clonotype)].groupby([self.clonotype_fmt]).binding_concordance.max()
        pep_concordance_scores = self.sub_df.groupby([self.clonotype_fmt]).binding_concordance.max()
        relevant_clonotypes    = (all_concordance_scores == pep_concordance_scores).to_dict() # boolean dictionary
        #print(sub_df)
        #print(relevant_clonotypes)
        SpecificityMatrix.remaining_indexes.update(self.sub_df.index[~self.sub_clonotype.map(relevant_clonotypes)]) # Indexes of non-top-concordance

        high_concordance_data = self.sub_df[self.sub_clonotype.map(relevant_clonotypes)].sort_values(by=['gems_per_specificity','binding_concordance'], ascending=[False,False])
        non_redundant_indexes = high_concordance_data.drop_duplicates(subset=[self.clonotype_fmt, 'peptide_HLA','gems_per_specificity','binding_concordance']).index
        return non_redundant_indexes # Indexes of top-concordance  #, na_action='ignore'


# # MAIN
# # Peptide per GEM
def peptides_per_gem(df, clonotype_fmt="ct",
					 bc_threshold=1, tcr_threshold=1,
					 exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
					 show=True, save_tuba=False, save_sund=False, save_report=False):

    print('Inside peptides_per_gem')

    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], point_scale=6, grid_axis='y')
    fig.var2col(df[clonotype_fmt].astype(int).astype(str))

    for mhc_barcode, sub_df in df.sort_values(by=['epitope_rank', clonotype_fmt]).groupby('peptide_HLA', sort=False):
        gems = sub_df.gem
        mhcs = sub_df.peptide_HLA.values
        umis = sub_df.umi_count_mhc
        clon = sub_df[clonotype_fmt].astype(int).astype(str)
        cols = fig.get_color(clon)

        scatter = fig.ax.scatter(gems, mhcs, s=umis*fig.point_scale, c=cols, edgecolors='face', alpha=0.5) #, alpha=0.3

    fig.add_labels_and_ticks(by='GEM')
    fig.add_filtering_legend()
    fig.add_label_legend()
    fig.add_size_legend()
    fig.finalize(show, save_tuba, save_sund)


def multiple_peptides_per_gem(df, clonotype_fmt="ct",
                              bc_threshold=1, tcr_threshold=1,
                              exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
                              show=True, save_tuba=False, save_sund=False, save_report=False):

    print('Inside multiple_peptides_per_gem')

    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], point_scale=6, grid_axis='y')
    
    df['tcr_category'] = pd.Categorical(df.tcr_category, categories=['unique chains', 'multiple chains', 'missing chain'], ordered=True)

    xmin, xmax = -0.5, 0

    for section, (tcr_multiplet_cat, group) in enumerate(df.groupby('tcr_category')):
        # Loop over all peptides to plot current annotations for each GEM
        for mhc_barcode, sub_df in group.groupby('peptide_HLA', sort=False):
            gems = sub_df.gem
            mhcs = sub_df.peptide_HLA.values
            umis = sub_df.umi_count_mhc
            xmax += len(gems.unique())

            scatter = fig.ax.scatter(gems, mhcs, s=umis*fig.point_scale, c='grey', marker='1')

            if section % 2 == 0:
                plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)
            xmin = xmax-0.5

        # Loop over all peptides to plot 'hidden' annotations for each GEM
        for mhc_barcode, sub_df in group.groupby('peptide_HLA', sort=False):
            gems = sub_df.gem.to_list()
            mhcs = sub_df.peptide_HLA_lst.to_list()
            umis = sub_df.umi_count_lst_mhc.to_list()

            X = [gems[i] for i, data in enumerate(mhcs) for j in range(len(data))]
            Y = [val for data in mhcs for val in data]
            S = [val*fig.point_scale for data in umis for val in data]

            scatter = fig.ax.scatter(X, Y, s=S, edgecolors='face', alpha=0.5)

        # https://www.google.com/search?client=safari&rls=en&q=what+is+the+top+coordinate+of+my+plt+plot&ie=UTF-8&oe=UTF-8
        # https://matplotlib.org/3.2.0/tutorials/advanced/transforms_tutorial.html
        # https://riptutorial.com/matplotlib/example/16030/coordinate-systems-and-text
        fig.ax.text(xmax, len(df.peptide_HLA.unique())*1.1, tcr_multiplet_cat, horizontalalignment='right')
    fig.add_labels_and_ticks(by='GEM')
    fig.add_filtering_legend()
    fig.add_size_legend()
    fig.finalize(show, save_tuba, save_sund)

####################################

#def multiple_peptides_per_gem_w_filtering(df, clonotype_fmt="ct",
#                                          bc_threshold=1, tcr_threshold=1,
#                                          exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
#                                          show=True, save_tuba=False, save_sund=False, save_report=False):
#
#    print('Inside multiple_peptides_per_gem_w_filtering')
#
#    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], point_scale=6, grid_axis='y')
#    
#    df['tcr_category'] = pd.Categorical(df.tcr_category, categories=['unique chains', 'multiple chains', 'missing chain'], ordered=True)
#
#    xmin, xmax = -0.5, 0
#
#    for section, (tcr_multiplet_cat, group) in enumerate(df.groupby('tcr_category')):
#        # Loop over all peptides to plot current annotations for each GEM
#        for mhc_barcode, sub_df in group.groupby('peptide_HLA', sort=False):
#            gems = sub_df.gem
#            mhcs = sub_df.peptide_HLA.values
#            umis = sub_df.umi_count_mhc
#            xmax += len(gems.unique())
#
#            scatter = fig.ax.scatter(gems, mhcs, s=umis*fig.point_scale, c='grey', marker='1')
#
#            if section % 2 == 0:
#                plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)
#            xmin = xmax-0.5
#
#        # Loop over all peptides to plot 'hidden' annotations for each GEM
#        for mhc_barcode, sub_df in group.groupby('peptide_HLA', sort=False):
#            gems = sub_df.gem.to_list()
#            mhcs = sub_df.peptide_HLA_lst.to_list()
#            umis = sub_df.umi_count_lst_mhc.to_list()
#
#            X = np.array([gems[i] for i, data in enumerate(mhcs) for j in range(len(data))])
#            Y = np.array([val for data in mhcs for val in data])
#            S = np.array([val*fig.point_scale for data in umis for val in data])
#
#            B = S >= bc_threshold * fig.point_scale
#
#            X = X[B]
#            Y = Y[B]
#            S = S[B]
#
#            if exclude_specificity_singlets:
#                idxs = get_nonsinglet_idxs(Y)
#                X = X[idxs]
#                Y = Y[idxs]
#                S = S[idxs]
#
#            scatter = fig.ax.scatter(X, Y, s=S, edgecolors='face', alpha=0.5)
#
#        # https://www.google.com/search?client=safari&rls=en&q=what+is+the+top+coordinate+of+my+plt+plot&ie=UTF-8&oe=UTF-8
#        # https://matplotlib.org/3.2.0/tutorials/advanced/transforms_tutorial.html
#        # https://riptutorial.com/matplotlib/example/16030/coordinate-systems-and-text
#        fig.ax.text(xmax, len(df.peptide_HLA.unique())*1.1, tcr_multiplet_cat, horizontalalignment='right')
#    fig.add_labels_and_ticks(by='GEM')
#    fig.add_filtering_legend()
#    fig.add_size_legend()
#    fig.finalize(show, save_tuba, save_sund)

####################################

def multiple_peptides_per_gem_w_filtering(df, clonotype_fmt="ct",
                                          bc_threshold=1, tcr_threshold=1,
                                          threshold=pd.Series(),
                                          exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
                                          show=True, save_tuba=False, save_sund=False, save_report=False):

    print('Inside multiple_peptides_per_gem_w_filtering')

    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], grid_axis='y') #point_scale=6
    
    df['tcr_category'] = pd.Categorical(df.tcr_category, categories=['unique chains', 'multiple chains', 'missing chain'], ordered=True)

    xmin, xmax = -0.5, 0

    for section, (tcr_multiplet_cat, group) in enumerate(df.groupby('tcr_category')):
        # Loop over all peptides to plot current annotations for each GEM
        for mhc_barcode, sub_df in group.groupby('peptide_HLA', sort=False):
            gems = sub_df.gem
            mhcs = sub_df.peptide_HLA.values
            umis = sub_df.umi_count_mhc
            xmax += len(gems.unique())

            scatter = fig.ax.scatter(gems, mhcs, s=umis*fig.point_scale, c='grey', marker='.', alpha=0) #c='grey', marker='1'

            if section % 2 == 0:
                plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)
            xmin = xmax-0.5

        # Loop over all peptides to plot 'hidden' annotations for each GEM
        for mhc_barcode, sub_df in group.groupby('peptide_HLA', sort=False):
            expl_df = (sub_df.drop(columns=['peptide_HLA','umi_count_mhc','umi_count_TRA','umi_count_TRB'])
                             .explode(['umi_count_lst_mhc','peptide_HLA_lst'])
                             .explode('umi_count_lst_TRA')
                             .explode('umi_count_lst_TRB')
                             .rename(columns={'peptide_HLA_lst':'peptide_HLA',
                                              'umi_count_lst_mhc':'umi_count_mhc',
                                              'umi_count_lst_TRA':'umi_count_TRA',
                                              'umi_count_lst_TRB':'umi_count_TRB'}))

            if not threshold.empty:
                filtering = eval(' & '.join(['(expl_df.%s >= %d)'%(k,v) for k,v in threshold.items()]))
                expl_df = expl_df[filtering].copy()

            gems = expl_df.gem
            mhcs = expl_df.peptide_HLA
            umis = expl_df.umi_count_mhc.astype(float)

            scatter = fig.ax.scatter(gems, mhcs, s=umis*fig.point_scale, edgecolors='face', alpha=0.5)

        # https://www.google.com/search?client=safari&rls=en&q=what+is+the+top+coordinate+of+my+plt+plot&ie=UTF-8&oe=UTF-8
        # https://matplotlib.org/3.2.0/tutorials/advanced/transforms_tutorial.html
        # https://riptutorial.com/matplotlib/example/16030/coordinate-systems-and-text
        fig.ax.text(xmax, len(df.peptide_HLA.unique())*1.1, tcr_multiplet_cat, horizontalalignment='right')
    fig.add_labels_and_ticks(by='GEM')
    fig.add_filtering_legend()
    #fig.add_size_legend()
    fig.finalize(show, save_tuba, save_sund)

####################################
####################################

def multiple_peptides_per_gem_w_filtering(df, clonotype_fmt="ct", y='peptide_HLA', color_by_hashing=None,
                                          bc_threshold=1, tcr_threshold=1,
                                          threshold=pd.Series(),
                                          exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
                                          show=True, save_tuba=False, save_sund=False, save_report=False):

    print('Inside multiple_peptides_per_gem_w_filtering')
    if y == 'peptide_HLA':
        by = 'GEM'
        y_lst = 'peptide_HLA_lst'
        umi = 'umi_count_mhc'
        umi_lst = 'umi_count_lst_mhc'
        sort_values = ['epitope_rank', 'HLA_match', umi] #clonotype_fmt
        scale = 6
    elif y == 'sample_id':
        by = 'somethingelse'
        y_lst = 'sample_id_lst'
        umi = 'umi_count_cd8'
        umi_lst = 'umi_count_lst_cd8'
        sort_values = [y, 'HLA_match', umi] #'sample_id'
        scale = 1

    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], grid_axis='y',point_scale=scale) #point_scale=6
    
    df['tcr_category'] = pd.Categorical(df.tcr_category, categories=['unique chains', 'multiple chains', 'missing chain'], ordered=True)

    xmin, xmax = -0.5, 0

    for section, (tcr_multiplet_cat, group) in enumerate(df.groupby('tcr_category')):
        # Loop over all peptides to plot current annotations for each GEM
        for mhc_barcode, sub_df in group.sort_values(by=sort_values).groupby(y, sort=False):
            gems = sub_df.gem
            mhcs = sub_df[y].values
            umis = sub_df[umi].values
            xmax += len(gems.unique())

            scatter = fig.ax.scatter(gems, mhcs, s=umis*fig.point_scale, c='grey', marker='.', alpha=0) #c='grey', marker='1'

            if section % 2 == 0:
                plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)
            xmin = xmax-0.5

        # Loop over all peptides to plot 'hidden' annotations for each GEM
        for mhc_barcode, sub_df in group.groupby(y, sort=False):
            expl_df = (sub_df.drop(columns=[y,umi]) #,'umi_count_TRA','umi_count_TRB'
                             .explode([umi_lst,y_lst])
                             #.explode('umi_count_lst_TRA')
                             #.explode('umi_count_lst_TRB')
                             .rename(columns={y_lst:y,
                                              umi_lst:umi})) #'umi_count_lst_TRA':'umi_count_TRA','umi_count_lst_TRB':'umi_count_TRB'
                                              
                                              

            if not threshold.empty:
                filtering = eval(' & '.join(['(expl_df.%s >= %d)'%(k,v) for k,v in threshold.items()]))
                expl_df = expl_df[filtering].copy()

            gems = expl_df.gem
            mhcs = expl_df[y].values
            umis = expl_df[umi].astype(float).values
            if color_by_hashing is not None:
                cols = expl_df.apply(lambda row: 'black' if row[y].split(' ')[-1] in row.HLA_cd8 else 'red', axis=1)
                scatter = fig.ax.scatter(gems, mhcs, s=umis*fig.point_scale, facecolors=cols, edgecolors='face', alpha=0.5)
            else:
                scatter = fig.ax.scatter(gems, mhcs, s=umis*fig.point_scale, edgecolors='face', alpha=0.5)

        # https://www.google.com/search?client=safari&rls=en&q=what+is+the+top+coordinate+of+my+plt+plot&ie=UTF-8&oe=UTF-8
        # https://matplotlib.org/3.2.0/tutorials/advanced/transforms_tutorial.html
        # https://riptutorial.com/matplotlib/example/16030/coordinate-systems-and-text
        fig.ax.text(xmax, len(df[y].unique())*1.1, tcr_multiplet_cat, horizontalalignment='right')
    fig.add_labels_and_ticks(by=by) #'GEM'
    fig.add_filtering_legend()
    fig.add_size_legend()
    fig.finalize(show, save_tuba, save_sund)

####################################

def peptide_per_clonotype_by_gem_size(df, clonotype_fmt="ct",
                                      bc_threshold=1, tcr_threshold=1,
                                      exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
                                      show=True, save_tuba=False, save_sund=False, save_report=False):

    fig = Plotting(bc_threshold, tcr_threshold, 'GEMs', df.gem, df[clonotype_fmt], point_scale=6)
  
    df = calc_binding_concordance(df.copy(), clonotype_fmt)

    for mhc_barcode, sub_df in df.sort_values(by=['epitope_rank']).groupby('peptide_HLA', sort=False):
        plateau = SpecificityMatrix(df, sub_df, mhc_barcode, clonotype_fmt)

        if clonotype_fmt == 'clonotype':
            tcrs = sub_df.loc[plateau.relevant_indexes, clonotype_fmt]
        else:
            tcrs = sub_df.loc[plateau.relevant_indexes, clonotype_fmt].astype(int).astype(str) #.str.split('clonotype', expand=True)[1]
        mhcs = sub_df.loc[plateau.relevant_indexes, 'peptide_HLA']
        gems = sub_df.loc[plateau.relevant_indexes, 'gems_per_specificity']
        cols = sub_df.loc[plateau.relevant_indexes, 'binding_concordance']
 
        scatter = fig.ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*fig.point_scale, label=None) #, viridis_r edgecolors='face', alpha=0.6, 
 
    print(SpecificityMatrix.remaining_indexes)
    tcrs = df.loc[SpecificityMatrix.remaining_indexes, clonotype_fmt].astype(int).astype(str)#.str.split('clonotype', expand=True)[1]
    mhcs = df.loc[SpecificityMatrix.remaining_indexes, 'peptide_HLA']
    gems = df.loc[SpecificityMatrix.remaining_indexes, 'gems_per_specificity']
    cols = df.loc[SpecificityMatrix.remaining_indexes, 'binding_concordance']
 
    scatter = fig.ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*fig.point_scale, label=None) #viridis_r edgecolors='face', edgecolors='None', alpha=0.6, 

    # Plotting decorations
    fig.add_labels_and_ticks()
    fig.add_colorbar()
    fig.add_filtering_legend()
    fig.add_size_legend()
    fig.finalize(show, save_tuba, save_sund)
 
    # Report lists of unique TCRs and GEMs
    if save_report:
        regex = re.findall('unique_tcrs', save_report[0])
        if regex:
            np.savetxt(save_report[0], fig.unique_tcrs, fmt='%s')
            np.savetxt(save_report[1], fig.unique_gems, fmt='%s')
        else:
            np.savetxt(save_report[1], fig.unique_tcrs, fmt='%s')
            np.savetxt(save_report[0], fig.unique_gems, fmt='%s')


#########################
####################################

def peptide_per_clonotype_by_gem_size(df, clonotype_fmt="ct", y='peptide_HLA', color_by_hashing=None,
                                      bc_threshold=1, tcr_threshold=1,
                                      exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
                                      show=True, save_tuba=False, save_sund=False, save_report=False):
    
    if y == 'peptide_HLA':
        by = 'GEM'
        y_lst = 'peptide_HLA_lst'
        umi = 'umi_count_mhc'
        umi_lst = 'umi_count_lst_mhc'
        sort_values = ['epitope_rank','binding_concordance']
    elif y == 'sample_id':
        by = 'somethingelse'
        y_lst = 'sample_id_lst'
        umi = 'umi_count_cd8'
        umi_lst = 'umi_count_lst_cd8'
        sort_values = ['sample_id', 'binding_concordance']

    fig = Plotting(bc_threshold, tcr_threshold, 'GEMs', df.gem, df[clonotype_fmt], point_scale=6)
    
    df = calc_binding_concordance(df.copy(), clonotype_fmt)
    
    if color_by_hashing:
        concordance = 'hla_concordance'
    else:
        concordance = 'binding_concordance'
        
    for mhc_barcode, sub_df in df.sort_values(by=sort_values).groupby(y, sort=False):
        plateau = SpecificityMatrix(df, sub_df, mhc_barcode, clonotype_fmt)

        if clonotype_fmt == 'clonotype':
            tcrs = sub_df.loc[plateau.relevant_indexes, clonotype_fmt]
        else:
            tcrs = sub_df.loc[plateau.relevant_indexes, clonotype_fmt].astype(int).astype(str) #.str.split('clonotype', expand=True)[1]
        mhcs = sub_df.loc[plateau.relevant_indexes, y]
        gems = sub_df.loc[plateau.relevant_indexes, 'gems_per_specificity']
        cols = sub_df.loc[plateau.relevant_indexes, concordance] #'binding_concordance'
 
        scatter = fig.ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*fig.point_scale, label=None) #, viridis_r edgecolors='face', alpha=0.6, 
 
    print(SpecificityMatrix.remaining_indexes)
    tcrs = df.loc[SpecificityMatrix.remaining_indexes, clonotype_fmt].astype(int).astype(str) #.str.split('clonotype', expand=True)[1] 
    mhcs = df.loc[SpecificityMatrix.remaining_indexes, y]
    gems = df.loc[SpecificityMatrix.remaining_indexes, 'gems_per_specificity']
    cols = df.loc[SpecificityMatrix.remaining_indexes, concordance] #'binding_concordance'
 
    scatter = fig.ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*fig.point_scale, label=None) #viridis_r edgecolors='face', edgecolors='None', alpha=0.6, 

    # Plotting decorations
    fig.add_labels_and_ticks(by='somethingelse') #OBS!
    fig.add_colorbar()
    fig.add_filtering_legend()
    fig.add_size_legend()
    fig.finalize(show, save_tuba, save_sund)
    SpecificityMatrix.remaining_indexes.clear()
 
    # Report lists of unique TCRs and GEMs
    if save_report:
        regex = re.findall('unique_tcrs', save_report[0])
        if regex:
            np.savetxt(save_report[0], fig.unique_tcrs, fmt='%s')
            np.savetxt(save_report[1], fig.unique_gems, fmt='%s')
        else:
            np.savetxt(save_report[1], fig.unique_tcrs, fmt='%s')
            np.savetxt(save_report[0], fig.unique_gems, fmt='%s')


#########################

def peptide_per_clonotype_by_umi_counts(df, clonotype_fmt="ct",
                                        bc_threshold=1, tcr_threshold=1,
                                        exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
                                        show=True, save_tuba=False, save_sund=False, save_report=False):

    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], point_scale=2)

    df = calc_binding_concordance(df.copy(), clonotype_fmt)

    for mhc_barcode, sub_df in df.sort_values(by=['epitope_rank']).groupby('peptide_HLA', sort=False):
        plateau = SpecificityMatrix(df, sub_df, mhc_barcode, clonotype_fmt)

        clonotype_grp = sub_df.loc[plateau.relevant_indexes].groupby(clonotype_fmt)

        gems = clonotype_grp.umi_count_mhc.mean()
        errs = clonotype_grp.umi_count_mhc.std()
        cols = clonotype_grp.binding_concordance.mean()
        tcrs = gems.index.get_level_values(0).astype(int).astype(str)
        mhcs = [mhc_barcode] * len(tcrs)

        fig.ax.scatter(tcrs, mhcs, c='red', s=gems+errs, edgecolors='face',alpha=0.3)
        scatter = fig.ax.scatter(tcrs, mhcs, facecolor=cols, edgecolors=None, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*fig.point_scale, label=None) #, edgecolors='face'

    clonotype_grp = df.loc[plateau.remaining_indexes].groupby([clonotype_fmt, 'peptide_HLA'])

    gems = clonotype_grp.umi_count_mhc.mean()
    errs = clonotype_grp.umi_count_mhc.std()
    cols = clonotype_grp.binding_concordance.mean()
    tcrs = gems.index.get_level_values(0).astype(int).astype(str)
    mhcs = gems.index.get_level_values(1)

    fig.ax.scatter(tcrs, mhcs, c='red', s=gems+errs, edgecolors='face',alpha=0.3)
    scatter = fig.ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*fig.point_scale, label=None) #, edgecolors='face'

    # Plotting decorations
    fig.add_labels_and_ticks()
    fig.add_colorbar()
    fig.add_filtering_legend()
    fig.add_size_legend()
    fig.finalize(show, save_tuba, save_sund)
 
    # Report lists of unique TCRs and GEMs

    if save_report:
        regex = re.findall('unique_tcrs', save_report[0])
        if regex:
            np.savetxt(save_report[0], fig.unique_tcrs, fmt='%s')
            np.savetxt(save_report[1], fig.unique_gems, fmt='%s')
        else:
            np.savetxt(save_report[1], fig.unique_tcrs, fmt='%s')
            np.savetxt(save_report[0], fig.unique_gems, fmt='%s')

#########################

def mhc_read_count_per_clonotype(df, clonotype_fmt="ct",
                                 bc_threshold=1, tcr_threshold=1,
                                 exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
                                 show=True, save_tuba=False, save_sund=False, save_report=False):

    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], point_scale=2, grid_axis='y')
    fig.var2col(df['peptide_HLA'])

    xmin, xmax = -0.5, 0

    for i, (clonotype, sub_df) in enumerate(df.sort_values(by=['umi_count_mhc', 'umi_count_tcr'], ascending=False).groupby(clonotype_fmt, sort=True)):
        for peptide_HLA, s_df in sub_df.groupby('peptide_HLA', sort=False):
            gems = s_df.gem
            umis = s_df.umi_count_mhc
            peps = s_df.peptide_HLA
            cols = fig.get_color(peps)

            scatter = fig.ax.scatter(gems, umis, c=cols)

            xmax += len(np.unique(gems))
        if i % 2 == 0:
            plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)
        xmin = xmax-0.5

    fig.add_labels_and_ticks(by='GEM', ylabel="MHC BC UMI count")
    if df.umi_count_mhc.max() == df.umi_count_mhc.max(): # Test if max is NaN
        plt.ylim(0, int(df.umi_count_mhc.max())+1)
    fig.add_filtering_legend()
    fig.add_label_legend(columns=7, xanchor=0.02)
    fig.finalize(show, save_tuba, save_sund)

#########################
# plot boxplots instead? https://stackoverflow.com/questions/29779079/adding-a-scatter-of-points-to-a-boxplot-using-matplotlib
def tcr_read_count_per_clonotype(df, clonotype_fmt="ct",
                                 bc_threshold=1, tcr_threshold=1,
                                 exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
                                 show=True, save_tuba=False, save_sund=False, save_report=False):
    
    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], point_scale=2, grid_axis='y')
    fig.var2col(df['peptide_HLA'])

    box = list()
    scx = list()
    lab = list()
    col = list()
    xmx = list()

    for i, (clonotype, sub_df) in enumerate(df.sort_values(by=['umi_count_tcr', 'umi_count_mhc']).groupby(clonotype_fmt, sort=True)):
        for peptide_HLA, s_df in sub_df.groupby('peptide_HLA', sort=False):
            gems = s_df.gem
            umis = s_df.umi_count_tcr
            peps = s_df.peptide_HLA
            cols = fig.get_color(peps)

            box.append(umis.to_list())
            scx.append(np.random.normal(0, 0.04, len(umis)))
            lab.append(peptide_HLA)
            col.append(cols[0])
            xmx.append(i)

            #scatter = fig.ax.scatter(gems, umis, c=cols)

    xmin, xmax = -0.5, 0

    if len(lab) > 0:
        bp = fig.ax.boxplot(box, labels=lab, showfliers=False)
        
    for x, val, c, i in zip(scx, box, col, xmx):
        xmax += 1
        fig.ax.scatter(x+xmax, val, c=c, alpha=0.8)

        if i % 2 == 0:
            plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)
        xmin = xmax-0.5

    fig.add_labels_and_ticks(by='alt', ylabel="TCR UMI count")
    if df.umi_count_tcr.max() == df.umi_count_tcr.max(): # Test if max is NaN
        plt.ylim(0, int(df.umi_count_tcr.max())+1)
    fig.add_filtering_legend()
    #fig.add_label_legend(columns=7, xanchor=0.02)
    fig.finalize(show, save_tuba, save_sund)

#########################

## plot boxplots instead? https://stackoverflow.com/questions/29779079/adding-a-scatter-of-points-to-a-boxplot-using-matplotlib
#def tcr_read_count_per_clonotype(df, clonotype_fmt="ct",
#                                 bc_threshold=1, tcr_threshold=1,
#                                 exclude_clonotype_singlets=False, exclude_specificity_singlets=False,
#                                 show=True, save_tuba=False, save_sund=False, save_report=False):
#    
#    fig = Plotting(bc_threshold, tcr_threshold, 'UMIs', df.gem, df[clonotype_fmt], point_scale=2, grid_axis='y')
#    fig.var2col(df['peptide_HLA'])
#
#    xmin, xmax = -0.5, 0
#
#    for i, (clonotype, sub_df) in enumerate(df.sort_values(by=['umi_count_tcr', 'umi_count_mhc']).groupby(clonotype_fmt, sort=True)):
#        gems = s_df.gem
#        umis = s_df.umi_count_tcr
#        peps = s_df.peptide_HLA
#        cols = fig.get_color(peps)
#
#        scatter = fig.ax.scatter(gems, umis, c=cols)
#
#        xmax += len(np.unique(gems))
#        if i % 2 == 0:
#            plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)
#        xmin = xmax-0.5
#
#    fig.add_labels_and_ticks(by='GEM', ylabel="TCR UMI count")
#    if df.umi_count_tcr.max() == df.umi_count_tcr.max(): # Test if max is NaN
#        plt.ylim(0, int(df.umi_count_tcr.max())+1)
#    fig.add_filtering_legend()
#    fig.add_label_legend(columns=7, xanchor=0.02)
#    fig.finalize(show, save_tuba, save_sund)


#########################
