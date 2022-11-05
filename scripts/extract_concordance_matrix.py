#!/usr/bin/env python
# coding: utf-8

# In[1]:

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import itertools

plt.style.use('ggplot')

# # Functions
def calc_binding_concordance_archived(df):
    assert df.size > 0, "df empty"
    gems_per_specificity_df = df.groupby(['clonotype','epitope']).gem.count().to_frame().reset_index()
    gems_per_specificity_df.rename(columns={'gem': 'gems_per_specificity'}, inplace=True)
    
    gems_per_clonotype_df = df.groupby(['clonotype']).gem.count().to_frame().reset_index()
    gems_per_clonotype_df.rename(columns={'gem': 'gems_per_clonotype'}, inplace=True)
    
    df = pd.merge(df, gems_per_specificity_df, on=['clonotype', 'epitope'], how='left', suffixes=('_total', '')).merge(gems_per_clonotype_df, on='clonotype', how='left', suffixes=('_total', ''))
    df['binding_concordance'] = df.gems_per_specificity / df.gems_per_clonotype
    
    return df

def calc_binding_concordance(df, clonotype_fmt):
    #assert df.size > 0, "df empty"
    gems_per_specificity = df.groupby([clonotype_fmt,'epitope']).gem.count().to_dict()
    df['gems_per_specificity'] = df.set_index([clonotype_fmt,'epitope']).index.map(gems_per_specificity)

    gems_per_clonotype = df.groupby([clonotype_fmt]).gem.count().to_dict()
    df['gems_per_clonotype'] = df[clonotype_fmt].map(gems_per_clonotype)
    
    df['binding_concordance'] = df.gems_per_specificity / df.gems_per_clonotype
    
    #return df

def epitope_sorter_index(df):
    EPITOPE_SORTER = ['CLYBL',
                      'v9', 'v15', 'v19', 'v3', 'v5', 'v6', 'v10', 'v13',
                      'v16', 'v24', 'v41', 'v2', 'v11', 'v18', 'v23', 'v25',
                      'v26', 'v27', 'v4', 'v7', 'v8', 'v12', 'v14', 'v20',
                      'v21', 'v1', 'v17', 'v22' ,'v30', 'v31', 'v36', 'v37',
                      'v38', 'v39', 'v40', 'v32', 'v33', 'v34', 'v35']

    sorterIndex = dict(zip(EPITOPE_SORTER,range(len(EPITOPE_SORTER))))
    
    return df.epitope.map(sorterIndex) #df['epitope_rank'] = 

def filter_levels(df):
    filter_options = [(df.copy(), 'All TCR annotations', 'no_filtration'),
                      (df[df.single_chain_only == False].copy(), 'TCRs with minimum 1 annotation for each chain', 'exclude_single-chain_TCRs'),
                      (df[(df.single_chain_only == False) &
                          (df.single_TRA == True) &
                          (df.single_TRB == True)].copy(), 'TCRs with unique annotations for each chain', 'exclude_ambiguous_and_single-chain_TCRs')]
    for plot_df, title, name in filter_options:
        yield plot_df, title, name

def filter_low_level(credible_df, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets, clonotype_fmt):
	flt_df = credible_df[(credible_df.umi_count_mhc >= read_threshold) &
                         (credible_df.umi_count_tcr >= tcr_threshold)].copy()
	print('UMI filtering', flt_df.gem.unique().shape[0])
	if exclude_clonotype_singlets:
		flt_df.drop(flt_df[~flt_df.duplicated(subset=clonotype_fmt, keep=False)].index, inplace=True)
	if exclude_specificity_singlets:
		flt_df.drop(flt_df[~flt_df.duplicated(subset=[clonotype_fmt, 'epitope'], keep=False)].index, inplace=True)
	return flt_df

def zip_tcr_concordance(tcrs, cols):
    tcr2concordance = dict()
    for tcr in tcrs.unique():
        indexes = np.where(tcrs.values == tcr)
        tcr2concordance[tcr] = cols.values[indexes]
    return tcr2concordance

# # MAIN
# # Peptide per GEM
def peptides_per_gem(credible_df, version=1, show=True, save_tuba=False, save_sund=False):
    """Make specificity plots"""
    
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'y'
    
    sortby = 'num_clonotype'
    credible_df.sort_values(by=['epitope_rank', sortby], inplace=True)
    
    # GEM to clonotype
    gem_to_clonotype = dict()
    for gem in credible_df.gem.unique():
        clonotypes = credible_df[credible_df.gem == gem].clonotype.values
        assert len(np.unique(clonotypes)) == 1, print(clonotypes, gem)
        gem_to_clonotype[gem] = clonotypes[0]

    # Clonotype to color
    all_clonotypes = credible_df.clonotype.unique()
    col_clonotypes = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'] * len(all_clonotypes)
    clonotype_to_color = dict()
    for i, clonotype in enumerate(all_clonotypes):
        clonotype_to_color[clonotype] = col_clonotypes[i]
    
    # Plotting
    project = "peptide_per_gem/"
    
    for read_threshold in [1, 2, 4, 10]: #, 20, 50
        for tcr_threshold in [1, 2, 4, 10]: #, 20, 50
            for exclude_clonotype_singlets in [False]: #, True
                for exclude_specificity_singlets in [False]: #, True
                    unique_gems = set()

                    fig, ax = plt.subplots(figsize=(20, 10))
                    for mhc_barcode in credible_df.peptide_HLA.unique():
                        sub_df = credible_df[(credible_df.peptide_HLA == mhc_barcode) &
                                             (credible_df.umi_count_mhc >= read_threshold) &
                                             (credible_df.umi_count_tcr >= tcr_threshold)].copy()
                        if exclude_clonotype_singlets:
                            sub_df.drop(sub_df[~sub_df.duplicated(subset='num_clonotype', keep=False)].index, inplace=True)
                        if exclude_specificity_singlets:
                            sub_df.drop(sub_df[~sub_df.duplicated(subset=['num_clonotype', 'epitope'], keep=False)].index, inplace=True)

                        gems = sub_df.gem.values
                        mhcs = [mhc_barcode] * len(gems)
                        umis = sub_df.umi_count_mhc.values

                        clonotypes = np.array([gem_to_clonotype[gem] for gem in gems])
                        colors = [clonotype_to_color[ct] for ct in clonotypes]
                        unique_gems.update(sub_df.gem.to_list())

                        scatter = ax.scatter(gems, mhcs, s=umis, c=colors, edgecolors='face') #, alpha=0.3

                    plt.xlabel("%i GEMs" %len([item.get_text() for item in ax.get_xticklabels()]), fontsize=16)
                    plt.tick_params(labelbottom=False, labelright=True, labelsize=8) #8

                    # Criteria
                    textstr = '\n'.join((
                        "Criteria",
                        "Min. BC UMI count \t %i".expandtabs() %read_threshold,
                        "Min. TCR UMI count \t\t    %i".expandtabs() %tcr_threshold,
                        "Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
                        "Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
                    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
                    #textstr = '\n'.join(("test \t test".expandtabs(),"test \t %i".expandtabs() %(read_threshold),"tester \t test".expandtabs()))
                    ax.text(0.05, 0.95, textstr,
                            fontsize=12,
                            horizontalalignment='left',
                            verticalalignment='top',
                            clip_on=False,
                            transform=ax.transAxes,
                            bbox=props)

                    from matplotlib.lines import Line2D

                    legend_elements = []
                    for clonotype_label in all_clonotypes:
                        legend_elements += [Line2D([0], [0], marker='o', color='w', label=clonotype_label, markerfacecolor=clonotype_to_color[clonotype_label], markersize=10)]

                    legend1 = ax.legend(handles=legend_elements, ncol=9, loc=2, bbox_to_anchor=(0.02, -0.03))
                    ax.add_artist(legend1)

                    for umi_size in [1, 5, 10]:
                        plt.scatter([], [], c='k', alpha=0.3, s=umi_size, label=str(umi_size) + ' UMIs')
                    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Barcode UMIs', loc='lower right')


                    plt.ylabel("Peptide, HLA", fontsize=16)
                    plt.title("Peptide specificity per GEM", fontsize=20)
                    
                    if show:
                        plt.show()
                        print("OBS! Figures at not saved!")
                        return
                    if save_tuba:
                    	# OBS! save_tuba + project + 
                        plt.savefig("v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(3, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    if save_sund:
                        plt.savefig(save_sund + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    
                    plt.cla()   # Clear axis
                    plt.clf()   # Clear figure
                    plt.close(fig) # Close a figure window

def initialize_plotting_frame(ax, x,y):
	sentinel, = ax.plot(x, np.append(y, [y[0]]*(len(x)-len(y))))
	sentinel.remove()

def multiple_peptides_per_gem(credible_df, clonotype_fmt="ct", filtration="no_filtration", show=True, save_tuba=False, save_sund=False, save_report=False):
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'y'

    credible_df.sort_values(by=['epitope_rank'], inplace=True)
    print('Sorted dataframe', credible_df.gem.unique().shape[0])

    project = "multiple_peptides_per_gem"

    xmin, xmax = -0.5, 0
    
    credible_df['tcr_category'] = pd.Categorical(credible_df.tcr_category, categories=['unique chains', 'missing chain', 'multiple chains'], ordered=True)

    fig, ax = plt.subplots(figsize=(20, 10))
    for section, tcr_multiplet_cat in enumerate(credible_df.tcr_category.cat.categories):
        for mhc_barcode in credible_df.peptide_HLA.unique():
            sub_df = credible_df[(credible_df.tcr_category == tcr_multiplet_cat) & (credible_df.peptide_HLA == mhc_barcode)]
            if sub_df.shape[0] == 0:
                continue

            gems = sub_df.gem.values
            mhcs = [mhc_barcode] * len (gems)
            umis = sub_df.umi_count_mhc.values * 5
            xmax += len(np.unique(gems))

            scatter = ax.scatter(gems, mhcs, s=umis, c='grey', marker='1')

            if section % 2 == 0:
                plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)

            #props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
            #ax.text((xmax-xmin)/2.0+xmin, 0.95, tcr_multiplet_cat,
            #        fontsize=12,
            #        horizontalalignment='left',
            #        verticalalignment='top',
            #        clip_on=False,
            #        transform=ax.transAxes,
            #        bbox=props)
            xmin = xmax-0.5

        for mhc_barcode in credible_df.peptide_HLA.unique():
            sub_df = credible_df[(credible_df.tcr_category == tcr_multiplet_cat) & (credible_df.peptide_HLA == mhc_barcode)]
            if sub_df.shape[0] == 0:
                continue

            gems = sub_df.gem.values
            mhcs = sub_df.peptide_HLA_lst.to_list()
            umis = sub_df.umi_count_lst_mhc.to_list()

            X = [gems[i] for i, data in enumerate(mhcs) for j in range(len(data))]
            Y = [val for data in mhcs for val in data]
            S = [val*5 for data in umis for val in data]

            scatter = ax.scatter(X, Y, s=S, edgecolors='face', alpha=0.5)

        # https://www.google.com/search?client=safari&rls=en&q=what+is+the+top+coordinate+of+my+plt+plot&ie=UTF-8&oe=UTF-8
        # https://matplotlib.org/3.2.0/tutorials/advanced/transforms_tutorial.html
        # https://riptutorial.com/matplotlib/example/16030/coordinate-systems-and-text
        ax.text(xmax, len(credible_df.peptide_HLA.unique())*1.1, tcr_multiplet_cat, horizontalalignment='right')

    number_gems = len([item.get_text() for item in ax.get_xticklabels()])
    plt.xlabel("%i GEMs" %number_gems, fontsize=16)
    plt.tick_params(labelbottom=False, labelright=True, labelsize=8) #8
    plt.xlim(-1, number_gems+1)
    for umi_size in [1, 5, 10]:
        plt.scatter([], [], c='k', alpha=0.3, s=umi_size*5, label=str(umi_size) + ' UMIs')
        plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Barcode UMIs', loc='upper left')
        plt.ylabel("Peptide, HLA", fontsize=16)
        plt.title("Peptide specificity per GEM sorted by TCR chain annotations", fontsize=20)

    if show:
        plt.show()
        print("OBS! Figures are not saved!")
        return
    if save_tuba:
        plt.savefig(project + '_' + save_tuba +'.png')



def peptide_per_clonotype_by_gem_size(credible_df, clonotype_fmt="ct", filtration="no_filtration", bc_threshold=1, tcr_threshold=1, exclude_clonotype_singlets=False, exclude_specificity_singlets=False, show=True, save_tuba=False, save_sund=False, save_report=False):
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'both'#y
    
    #sortby = 'clonotype'
    credible_df.sort_values(by=['epitope_rank'], inplace=True) #,sortby
    print('Sorted dataframe', credible_df.gem.unique().shape[0])

    #filtrations = ['no_filtration', 'exclude_ambiguous_TCR_ann', 'exclude_single-chain_clonotypes', 'exclude_ambiguous_BC_ann', 'exclude_ambiguous_BC-TCR_and_single-chain_clonotypes']

    
    project = "peptide_per_clonotype_by_gem_size/"

    # OBS! Concordance 3D plots!
    #bc_umis = [1, 2, 3, 4, 5]
    #tcr_umis = [1, 2, 4, 6, 8, 10]
    #X, Y = np.meshgrid(bc_umis, tcr_umis)
    #z_concordance_matrix = np.zeros((len(tcr_umis), len(bc_umis)))
    #z_gem_matrix = np.zeros((len(tcr_umis), len(bc_umis)))
    
    #for x, bc_threshold in enumerate(bc_umis): # , 20, 50
    #    for y, tcr_threshold in enumerate(tcr_umis): #, 20, 50
    #        for exclude_clonotype_singlets in [False]:
    #            for exclude_specificity_singlets in [False]:
                    # Subset the dataframe by the settings in the above for-loops
    flt_df = filter_low_level(credible_df, bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets, clonotype_fmt)
    if len(flt_df) == 0:
        print(clonotype_fmt, filtration, bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets)
        continue
    calc_binding_concordance(flt_df, clonotype_fmt)
    print('After filtering', flt_df.gem.unique().shape[0])
 
    # Initialize
    unique_gems, unique_tcrs, concordances_per_peptide, remaining_indexes = set(), list(), list(), set()
 
    fig, ax = plt.subplots(figsize=(20, 10))
    for mhc_barcode in credible_df.peptide_HLA.unique():
        sub_df = flt_df[flt_df.peptide_HLA == mhc_barcode]
 
        # Specific to this plotting function
        related_tcrs = sub_df[clonotype_fmt].values
        if len(related_tcrs) == 0:
        	continue
 
        # Select TCRs with maximum specificity concordance
        all_concordance_scores = flt_df[flt_df[clonotype_fmt].isin(related_tcrs)].groupby([clonotype_fmt]).binding_concordance.max()
        pep_concordance_scores = sub_df.groupby([clonotype_fmt]).binding_concordance.max()
        #sub_df = sub_df.loc[sub_df.clonotype.map((all_concordance_scores == pep_concordance_scores).to_dict())]
 
        relevant_clonotypes = (all_concordance_scores == pep_concordance_scores).to_dict() # boolean dictionary
        relevant_indexes = sub_df.index[sub_df[clonotype_fmt].map(relevant_clonotypes)] # Indexes of top-concordance  #, na_action='ignore'
        remaining_indexes.update(sub_df.index[~sub_df[clonotype_fmt].map(relevant_clonotypes)]) # Indexes of non-top-concordance
 
        tcrs = sub_df.loc[relevant_indexes, clonotype_fmt].astype(int).astype(str)#.str.split('clonotype', expand=True)[1]
        mhcs = sub_df.loc[relevant_indexes, 'peptide_HLA']
        gems = sub_df.loc[relevant_indexes, 'gems_per_specificity']
        cols = sub_df.loc[relevant_indexes, 'binding_concordance']
 
        scatter = ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*5, label=None) #, edgecolors='face'
 
        unique_gems.update(sub_df.loc[relevant_indexes, 'gem'])
        unique_tcrs.append(tcrs.drop_duplicates().to_list())
 
        # Combine TCRs and concordance measures into a dict
        concordances_per_peptide.append(zip_tcr_concordance(tcrs, cols))
 
    tcrs = flt_df.loc[remaining_indexes, clonotype_fmt].astype(int).astype(str)#.str.split('clonotype', expand=True)[1]
    mhcs = flt_df.loc[remaining_indexes, 'peptide_HLA']
    gems = flt_df.loc[remaining_indexes, 'gems_per_specificity']
    cols = flt_df.loc[remaining_indexes, 'binding_concordance']
 
    scatter = ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*5, alpha=0.5, label=None) #edgecolors='face', edgecolors='None', 
 
    assert len(unique_gems.intersection(set(flt_df.loc[remaining_indexes, 'gem']))) == 0
    unique_gems.update(flt_df.loc[remaining_indexes, 'gem'])
    unique_tcrs.append(tcrs.drop_duplicates().to_list())
    unique_tcrs = list(dict.fromkeys(list(itertools.chain.from_iterable(unique_tcrs)))) # Flatten list and remove duplicates while maintaining order
    assert len(unique_tcrs) == len(set(unique_tcrs))
 
    # Combine TCRs and concordance measures into a dict
    concordances_per_peptide.append(zip_tcr_concordance(tcrs, cols))
    # Calculate the sum of the maximum concordances for each clonotype divided by the number of clonotypes
    #z_concordance_matrix[y,x] = sum([max([value for concordance in concordances_per_peptide for value in concordance.get(key, [])]) for key in unique_tcrs]) / len(unique_tcrs)
    #z_gem_matrix[y,x] = len(unique_gems)
 
    plt.tick_params(labelsize=8) 
    plt.xticks(rotation=90, size=2)
 
    sm = plt.cm.ScalarMappable(cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1))
    if mpl.__version__ < '3.1':
        sm._A = []
    plt.colorbar(sm)
 
    # Criteria
    textstr = '\n'.join((
        "Criteria",
        "Min. BC UMI count\t%i".expandtabs() %bc_threshold,
        "Min. TCR UMI count\t%i".expandtabs() %tcr_threshold))#,
        #"Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
        #"Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
    ax.text(0.05, 0.95, textstr,
            fontsize=12,
            horizontalalignment='left',
            verticalalignment='top',
            clip_on=False,
            transform=ax.transAxes,
            bbox=props)
 
    for number_gems in [2, 10, 50]:
        plt.scatter([], [], c='k', alpha=0.3, s=number_gems*5, label=str(number_gems) + ' GEMs')
    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='GEMs', loc='lower right')
 
    plt.xlabel("%i Clonotypes (of %i GEMs)" %(len(unique_tcrs), len(unique_gems)), fontsize=16)
    plt.ylabel("Peptide, HLA", fontsize=16)
    plt.title("Specificity concordance per clonotype", fontsize=20)
 
    if show:
        plt.savefig("legend.pdf", bbox_inches='tight')
        plt.show()
        print("OBS! Figures are not saved!")
        return #unique_gems, unique_tcrs
    if save_tuba:
        #print(save_tuba + project + clonotype_fmt + "/%s/"%filtration + "b%i.t%i.ecs_%s.ess_%s.pdf" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets))
        #plt.savefig(save_tuba + project + clonotype_fmt + "/%s/"%filtration + "b%i.t%i.ecs_%s.ess_%s.pdf" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
        #print(save_tuba + "b%i.t%i.ecs_%s.ess_%s.pdf" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets))
        #plt.savefig(save_tuba + "b%i.t%i.ecs_%s.ess_%s.pdf" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
        plt.savefig(save_tuba, bbox_inches='tight')
    if save_sund:
        plt.savefig(save_sund + project + clonotype_fmt + "/%s/"%filtration + "b%i.t%i.ecs_%s.ess_%s.pdf" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close(fig) # Close a figure window
 
    # Report lists of unique TCRs and GEMs
    if save_report:
        #np.savetxt(save_report + project + clonotype_fmt + "/%s/"%filtration + "unique_tcrs/" + "b%i.t%i.ecs_%s.ess_%s.lst" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), unique_tcrs, fmt='%s')
        #np.savetxt(save_report + project + clonotype_fmt + "/%s/"%filtration + "unique_gems/" + "b%i.t%i.ecs_%s.ess_%s.lst" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), list(unique_gems), fmt='%s')
        regex = re.findall('unique_tcrs', save_report[0])
        if regex:
            np.savetxt(save_report[0], unique_tcrs, fmt='%s')
            np.savetxt(save_report[1], list(unique_gems), fmt='%s')
        else:
            np.savetxt(save_report[1], unique_tcrs, fmt='%s')
            np.savetxt(save_report[0], list(unique_gems), fmt='%s')
            
        #np.savetxt(save_report + "unique_tcrs/" + "b%i.t%i.ecs_%s.ess_%s.lst" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), unique_tcrs, fmt='%s')
        #np.savetxt(save_report + "unique_gems/" + "b%i.t%i.ecs_%s.ess_%s.lst" %(bc_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), list(unique_gems), fmt='%s')
    #if save_report:
    #    np.savetxt(save_report + "concordance_comparison_matrix/" + "X", X, fmt='%d')
    #    np.savetxt(save_report + "concordance_comparison_matrix/" + "Y", Y, fmt='%d')
    #    np.savetxt(save_report + "concordance_comparison_matrix/" + "Z_conc", z_concordance_matrix)
    #    np.savetxt(save_report + "concordance_comparison_matrix/" + "Z_gems", z_gem_matrix)

#########################

def peptide_per_clonotype_by_gem_size_archive(credible_df, version=1, show=True, save_tuba=False, save_sund=False, save_report=False):
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'both'#y
    
    sortby = 'clonotype'
    credible_df.sort_values(by=['epitope_rank', sortby], inplace=True)
    print(credible_df.gem.unique().shape[0])

    filtrations = ['no_filtration', 'exclude_ambiguous_TCR_ann', 'exclude_single-chain_clonotypes', 'exclude_ambiguous_BC_ann', 'exclude_ambiguous_BC-TCR_and_single-chain_clonotypes']

    
    project = "peptide_per_clonotype_by_gem_size/"
    
    for read_threshold in [1, 10]: # , 20, 50
        for tcr_threshold in [1, 10]: #, 20, 50
            for exclude_clonotype_singlets in [False, True]:
                for exclude_specificity_singlets in [False, True]:
                    # Subset the dataframe by the settings in the above for-loops
                    tmp_df = credible_df[(credible_df.umi_count_mhc >= read_threshold) &
                                         (credible_df.umi_count_tcr >= tcr_threshold)].copy()
                    print(tmp_df.gem.unique().shape[0])
                    if exclude_clonotype_singlets:
                        tmp_df.drop(tmp_df[~tmp_df.duplicated(subset='num_clonotype', keep=False)].index, inplace=True)
                    if exclude_specificity_singlets:
                        tmp_df.drop(tmp_df[~tmp_df.duplicated(subset=['num_clonotype', 'epitope'], keep=False)].index, inplace=True)
                    tmp_df = calc_binding_concordance_archived(tmp_df)
                    print(tmp_df.gem.unique().shape[0])

                    # Initialize
                    unique_gems, unique_tcrs, remaining_indexes = set(), list(), set()

                    fig, ax = plt.subplots(figsize=(20, 10))
                    for mhc_barcode in credible_df.peptide_HLA.unique():
                        sub_df = tmp_df[tmp_df.peptide_HLA == mhc_barcode]

                        related_tcrs = sub_df.clonotype.values
                        if len(related_tcrs) == 0:
                        	continue

                        # Select TCRs with maximum specificity concordance
                        all_concordance_scores = tmp_df[tmp_df.clonotype.isin(related_tcrs)].groupby(['clonotype']).binding_concordance.max()
                        pep_concordance_scores = sub_df.groupby(['clonotype']).binding_concordance.max()
                        #sub_df = sub_df.loc[sub_df.clonotype.map((all_concordance_scores == pep_concordance_scores).to_dict())]

                        relevant_clonotypes = (all_concordance_scores == pep_concordance_scores).to_dict() # boolean dictionary
                        relevant_indexes = sub_df.index[sub_df.clonotype.map(relevant_clonotypes)] # Indexes of top-concordance 
                        remaining_indexes.update(sub_df.index[~sub_df.clonotype.map(relevant_clonotypes)]) # Indexes of non-top-concordance

                        tcrs = sub_df.loc[relevant_indexes, 'num_clonotype'].astype(int).astype(str)#.str.split('clonotype', expand=True)[1]
                        mhcs = sub_df.loc[relevant_indexes, 'peptide_HLA']
                        gems = sub_df.loc[relevant_indexes, 'gems_per_specificity']
                        cols = sub_df.loc[relevant_indexes, 'binding_concordance']

                        scatter = ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*5, edgecolors='face')

                        unique_gems.update(sub_df.loc[relevant_indexes, 'gem'])
                        unique_tcrs.append(tcrs.drop_duplicates().to_list())

                    tcrs = tmp_df.loc[remaining_indexes, 'num_clonotype'].astype(int).astype(str)#.str.split('clonotype', expand=True)[1]
                    mhcs = tmp_df.loc[remaining_indexes, 'peptide_HLA']
                    gems = tmp_df.loc[remaining_indexes, 'gems_per_specificity']
                    cols = tmp_df.loc[remaining_indexes, 'binding_concordance']

                    scatter = ax.scatter(tcrs, mhcs, c=cols, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems*5, edgecolors='face')

                    assert len(unique_gems.intersection(set(tmp_df.loc[remaining_indexes, 'gem']))) == 0
                    unique_gems.update(tmp_df.loc[remaining_indexes, 'gem'])
                    unique_tcrs.append(tcrs.drop_duplicates().to_list())
                    unique_tcrs = list(dict.fromkeys(list(itertools.chain.from_iterable(unique_tcrs)))) # Flatten list and remove duplicates while maintaining order
                    assert len(unique_tcrs) == len(set(unique_tcrs))

                    plt.tick_params(labelsize=8) 
                    plt.xticks(rotation=90, size=2)

                    sm = plt.cm.ScalarMappable(cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1))
                    if mpl.__version__ < '3.1':
                        sm._A = []
                    plt.colorbar(sm)

                    # Criteria
                    textstr = '\n'.join((
                        "Criteria",
                        "Min. BC UMI count\t%i".expandtabs() %read_threshold,
                        "Min. TCR UMI count\t%i".expandtabs() %tcr_threshold))#,
                        #"Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
                        #"Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
                    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
                    ax.text(0.05, 0.95, textstr,
                            fontsize=12,
                            horizontalalignment='left',
                            verticalalignment='top',
                            clip_on=False,
                            transform=ax.transAxes,
                            bbox=props)

                    for number_gems in [2, 10, 50]:
                        plt.scatter([], [], c='k', alpha=0.3, s=number_gems*5, label=str(number_gems) + ' GEMs')
                    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='GEMs', loc='lower right')
                    

                    plt.xlabel("%i Clonotypes (of %i GEMs)" %(len(unique_tcrs), len(unique_gems)), fontsize=16)
                    plt.ylabel("Peptide, HLA", fontsize=16)
                    plt.title("Specificity concordance per clonotype", fontsize=20)

                    if show:
                        plt.show()
                        print("OBS! Figures are not saved!")
                        return unique_gems, unique_tcrs
                    if save_tuba:
                        plt.savefig(save_tuba + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    if save_sund:
                        plt.savefig(save_sund + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    plt.cla()   # Clear axis
                    plt.clf()   # Clear figure
                    plt.close(fig) # Close a figure window

                    # Report lists of unique TCRs and GEMs
                    if save_report:
	                    np.savetxt(save_report + project + "unique_tcrs/" + "v%i.b%i.t%i.ecs_%s.ess_%s.lst" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), unique_tcrs, fmt='%s')
	                    np.savetxt(save_report + project + "unique_gems/" + "v%i.b%i.t%i.ecs_%s.ess_%s.lst" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), list(unique_gems), fmt='%s')


# # Peptide per clonotype (GEM counts)

# In[6]:


def peptide_per_clonotype_by_gem_size_archive(credible_df, version=1, show=True, save_tuba=False, save_sund=False):
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'both'#y
    
    sortby = 'clonotype'
    credible_df.sort_values(by=['epitope_rank', sortby], inplace=True)
    
    project = "peptide_per_clonotype_by_gem_size/"
    
    for read_threshold in [1, 10]: # , 20, 50
        for tcr_threshold in [1, 10]: #, 20, 50
            for exclude_clonotype_singlets in [False, True]:
                for exclude_specificity_singlets in [False, True]:

                    unique_gems, unique_tcrs = set(), set()

                    fig, ax = plt.subplots(figsize=(20, 10))

                    for mhc_barcode in credible_df.peptide_HLA.unique():

                        sub_df = credible_df[(credible_df.read_count_mhc >= read_threshold) &
                                             (credible_df.umi_count_tcr >= tcr_threshold)].copy()
                        if exclude_clonotype_singlets:
                            sub_df.drop(sub_df[~sub_df.duplicated(subset='num_clonotype', keep=False)].index, inplace=True)
                        if exclude_specificity_singlets:
                            sub_df.drop(sub_df[~sub_df.duplicated(subset=['num_clonotype', 'epitope'], keep=False)].index, inplace=True)

                        sub_df = calc_binding_concordance_archived(sub_df)
                        sub_df = sub_df[sub_df.peptide_HLA == mhc_barcode]

                        tcrs = sub_df.clonotype.str.split('clonotype').str[1].unique()
                        mhcs = [mhc_barcode] * len(tcrs)
                        gems = sub_df.groupby(['clonotype']).gems_per_specificity.mean().values * 5
                        colors = sub_df.groupby(['clonotype']).binding_concordance.mean().values

                        unique_gems.update(sub_df.gem.to_list())
                        unique_tcrs.update(sub_df.clonotype.to_list())

                        scatter = ax.scatter(tcrs, mhcs, c=colors, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems, edgecolors='face')

                    plt.tick_params(labelsize=8) 
                    plt.xticks(rotation=90, size=2)

                    sm = plt.cm.ScalarMappable(cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1))
                    if mpl.__version__ < '3.1':
                        sm._A = []
                    plt.colorbar(sm)

                    # Criteria
                    textstr = '\n'.join((
                        "Criteria",
                        "Min. barcode read count \t %i".expandtabs() %read_threshold,
                        "Min. TCR read count \t\t    %i".expandtabs() %tcr_threshold,
                        "Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
                        "Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
                    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
                    ax.text(0.05, 0.95, textstr,
                            fontsize=12,
                            horizontalalignment='left',
                            verticalalignment='top',
                            clip_on=False,
                            transform=ax.transAxes,
                            bbox=props)

                    for number_gems in [2, 10, 50]:
                        plt.scatter([], [], c='k', alpha=0.3, s=number_gems*5, label=str(number_gems) + ' GEMs')
                    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='GEMs', loc='lower right')

                    plt.xlabel("%i Clonotypes (of %i GEMs)" %(len(unique_tcrs), len(unique_gems)), fontsize=16)
                    plt.ylabel("Peptide, HLA", fontsize=16)
                    plt.title("Specificity concordance per clonotype", fontsize=20)

                    if show:
                        plt.show()
                        print("OBS! Figures are not saved!")
                        return
                    if save_tuba:
                        plt.savefig(save_tuba + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    if save_sund:
                        plt.savefig(save_sund + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    plt.cla()   # Clear axis
                    plt.clf()   # Clear figure
                    plt.close(fig) # Close a figure window


# # Peptide per clonotype (read counts)

# In[7]:


def peptide_per_clonotype_read_counts(credible_df, version=1, show=True, save_tuba=False, save_sund=False):
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'both'#y
    
    #credible_df.epitope = credible_df.epitope.astype("category")
    #credible_df.epitope.cat.set_categories(EPITOPE_SORTER, inplace=True)
    sortby = 'clonotype'
    credible_df.sort_values(by=['epitope_rank', sortby], inplace=True)
    
    project = "peptide_per_clonotype_by_read_size/"
    
    exclude_low_concordance_clonotypes = False
    show_tcr_multiplicity = False
    
    for read_threshold in [1, 10, 20]: #, 20, 50
        for tcr_threshold in [1, 10, 20]: #, 20, 50
            for exclude_clonotype_singlets in [False, True]:
                for exclude_specificity_singlets in [False, True]:
                    #print(str(read_threshold), str(tcr_threshold), str(exclude_clonotype_singlets), str(exclude_specificity_singlets))

                    unique_gems, unique_tcrs = set(), set()

                    fig, ax = plt.subplots(figsize=(20, 10))

                    for mhc_barcode in credible_df.peptide_HLA.unique(): #credible_df.epitope.unique(): #specificity_df.columns:
                        sub_df = credible_df[(credible_df.read_count_mhc >= read_threshold) &
                                             (credible_df.umi_count_tcr >= tcr_threshold)].copy()
                        if exclude_clonotype_singlets:
                            sub_df.drop(sub_df[~sub_df.duplicated(subset='num_clonotype', keep=False)].index, inplace=True)
                        if exclude_specificity_singlets:
                            sub_df.drop(sub_df[~sub_df.duplicated(subset=['num_clonotype', 'epitope'], keep=False)].index, inplace=True)

                        sub_df = calc_binding_concordance_archived(sub_df)
                        sub_df = sub_df[sub_df.peptide_HLA == mhc_barcode]

                        if exclude_low_concordance_clonotypes:
                            sub_df.drop()

                        tcrs = sub_df.clonotype.str.split('clonotype').str[1].unique()
                        mhcs = [mhc_barcode] * len(tcrs)
                        gems = sub_df.groupby(['clonotype']).read_count_mhc.mean().values #groupby(['clonotype']).gem.count().values
                        errs = sub_df.groupby(['clonotype']).read_count_mhc.std().values
                        colors = sub_df.groupby(['clonotype']).binding_concordance.mean().values # should binding concordance be calculated after subsetting by the threshold?

                        unique_gems.update(sub_df.gem.to_list())
                        unique_tcrs.update(sub_df.clonotype.to_list())
                            
                        # How to show number of GEMs? Plot a different symbol if only one GEM?
                        ax.scatter(tcrs, mhcs, c='red', s=gems+errs, edgecolors='face',alpha=0.3) #, alpha=0.3
                        scatter = ax.scatter(tcrs, mhcs, c=colors, cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems, edgecolors='face') #, alpha=0.3

                        if show_tcr_multiplicity:
                            tcrs = sub_df[sub_df.single_tcell == False].clonotype.str.split('clonotype').str[1].unique()
                            mhcs = [mhc_barcode] * len(tcrs)
                            colors = np.where(sub_df[sub_df.single_tcell == False].groupby(['clonotype']).binding_concordance.mean().values > 0.5, 'white', 'k')
                            ax.scatter(tcrs, mhcs, c=colors, marker='+', edgecolors='k') #, alpha=0.3

                    plt.tick_params(labelsize=8) #labelbottom=False, labelright=True, 
                    plt.xticks(rotation=90, size=2)

                    #plt.colorbar(cmap='viridis_r')
                    sm = plt.cm.ScalarMappable(cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1))
                    if mpl.__version__ < '3.1':
                        sm._A = []
                    plt.colorbar(sm)

                    for number_gems in [10, 50, 100]:
                        plt.scatter([], [], c='k', alpha=0.3, s=number_gems, label=str(number_gems) + ' reads')
                    for marker, response in [('+', "TCR singlet")]: #('o', True), 
                        plt.scatter([], [], c='k', alpha=0.3, s=50, label=str(response), marker=marker)
                    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, loc='lower right') #title='Reads', 

                    # Criteria
                    # Criteria
                    textstr = '\n'.join((
                        "Criteria",
                        "Min. barcode read count \t %i".expandtabs() %read_threshold,
                        "Min. TCR read count \t\t    %i".expandtabs() %tcr_threshold,
                        "Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
                        "Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
                    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
                    ax.text(0.63, 0.146, textstr,
                            fontsize=12,
                            horizontalalignment='left',
                            verticalalignment='top',
                            clip_on=False,
                            transform=ax.transAxes,
                            bbox=props)

                    plt.xlabel("%i Clonotypes (of %i GEMs)" %(len(unique_tcrs), len(unique_gems)), fontsize=16)
                    plt.ylabel("Peptide, HLA", fontsize=16)
                    plt.title("Specificity concordance per clonotype", fontsize=20)
                    
                    if show:
                        plt.show()
                        print("OBS! Figures are not saved!")
                        return
                    if save_tuba:
                        plt.savefig(save_tuba + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    if save_sund:
                        plt.savefig(save_sund + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    plt.cla()   # Clear axis
                    plt.clf()   # Clear figure
                    plt.close(fig) # Close a figure window


# # MHC barcode read count per clonotype

# In[8]:


def mhc_read_count_per_clonotype(credible_df, version=1, show=True, save_tuba=False, save_sund=False):
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'y'#both

    # Epitope to color
    all_epitopes = credible_df.peptide_HLA.unique()
    col_epitopes = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'] * len(all_epitopes)
    epitope_to_color = dict()
    for i, epitope in enumerate(all_epitopes):
        epitope_to_color[epitope] = col_epitopes[i]
        
    # Detected response
    sortby = 'umi_count_mhc' #'umi_count_tcr' #
    credible_df.sort_values(by=['num_clonotype', sortby], inplace=True)
    
    project = "read_count_per_clonotype/"
    
    for read_threshold in [1, 2, 4, 10]: #, 20, 50
        for tcr_threshold in [1, 2, 4, 10]: #, 20, 50
            for exclude_clonotype_singlets in [True]: #False, 
                for exclude_specificity_singlets in [False]: #, True
                    unique_gems, unique_tcrs = set(), set()

                    fig, ax = plt.subplots(figsize=(20, 10))

                    xmin, xmax = -0.5, 0

                    for i, clonotype in enumerate(credible_df.num_clonotype.unique()):
                        sub_df = credible_df[(credible_df.num_clonotype == clonotype) &
                                             (credible_df.umi_count_mhc >= read_threshold) &
                                             (credible_df.umi_count_tcr >= tcr_threshold)].copy() # & (credible_df.clonotype != 'None') & (credible_df.epitope != '0')
                        if exclude_clonotype_singlets:
                            sub_df.drop(sub_df[~sub_df.duplicated(subset='num_clonotype', keep=False)].index, inplace=True)
                        if exclude_specificity_singlets:
                            sub_df.drop(sub_df[~sub_df.duplicated(subset=['num_clonotype', 'epitope'], keep=False)].index, inplace=True)

                        gems = sub_df.gem.to_list()
                        mhc_umi_count = sub_df.umi_count_mhc.to_list()

                        xmax += len(np.unique(gems))
                        unique_gems.update(sub_df.gem.to_list())
                        unique_tcrs.update(sub_df.clonotype.to_list())

                        epitopes = sub_df.peptide_HLA.to_list()
                        colors = [epitope_to_color[ep] for ep in epitopes]

                        # How to show number of GEMs? Plot a different symbol if only one GEM?
                        #ax.scatter(gems, [-5]*len(gems))
                        scatter = ax.scatter(gems, mhc_umi_count, c=colors) #, edgecolors='face', cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems, edgecolors='face', alpha=0.3

                        if i % 2 == 0:
                            plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)

                        xmin = xmax-0.5

                    from matplotlib.lines import Line2D

                    legend_elements = []
                    for epitope_label in all_epitopes:
                        legend_elements += [Line2D([0], [0], marker='o', color='w', label=epitope_label, markerfacecolor=epitope_to_color[epitope_label], markersize=10)]

                    legend1 = ax.legend(handles=legend_elements, ncol=7, loc=2, bbox_to_anchor=(0.02, -0.03))
                    ax.add_artist(legend1)

                    plt.tick_params(labelbottom=False, labelright=True, labelsize=8) #labelbottom=False, 
                    plt.xticks(rotation=90, size=2)

                    # Criteria
                    textstr = '\n'.join((
                        "Criteria",
                        "Min. BC UMI count \t %i".expandtabs() %read_threshold,
                        "Min. TCR UMI count \t\t    %i".expandtabs() %tcr_threshold,
                        "Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
                        "Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
                    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
                    ax.text(0.12, 0.98, textstr,
                            fontsize=12,
                            horizontalalignment='left',
                            verticalalignment='top',
                            clip_on=False,
                            transform=ax.transAxes,
                            bbox=props)

                    plt.xlabel("%i GEMs (sectioned per clonotype (%i))" %(len(unique_gems), len(unique_tcrs)), fontsize=16)
                    plt.ylabel("pMHC barcode UMI counts", fontsize=16)
                    plt.title("MHC barcode UMI counts per GEM per clonotype", fontsize=20)
                    
                    if show:
                        plt.show()
                        print("OBS! Figures at not saved!")
                        return
                    if save_tuba:
                    	# OBS! save_tuba + project + 
                        plt.savefig("v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(5, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    if save_sund:
                        plt.savefig(save_sund + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    plt.cla()   # Clear axis
                    plt.clf()   # Clear figure
                    plt.close(fig) # Close a figure window


# ## Response

# In[9]:


def mhc_read_count_per_clonotype_response(credible_df, version=1, show=True, save_tuba=False, save_sund=False):
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'y'#both

    # Epitope to color
    all_epitopes = credible_df.peptide_HLA.unique()
    col_epitopes = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'] * len(all_epitopes)
    epitope_to_color = dict()
    for i, epitope in enumerate(all_epitopes):
        epitope_to_color[epitope] = col_epitopes[i]
        
    # Detected response
    sortby = 'read_count_mhc' #'umi_count_tcr' #
    credible_df.sort_values(by=['num_clonotype', sortby], inplace=True)
    
    project = "read_count_per_clonotype_response/"
    
    for read_threshold in [1, 10]: #, 20, 50
        for tcr_threshold in [1, 10]: #, 20, 50
            for exclude_clonotype_singlets in [False, True]:
                for exclude_specificity_singlets in [False, True]:
                    unique_gems, unique_tcrs = set(), set()

                    fig, ax = plt.subplots(figsize=(20, 10))

                    xmin, xmax = -0.5, 0

                    for i, clonotype in enumerate(credible_df.num_clonotype.unique()):
                        for marker, response in [('o', True), ('+', False)]:
                            sub_df = credible_df[(credible_df.num_clonotype == clonotype) &
                                                 (credible_df.read_count_mhc >= read_threshold) &
                                                 (credible_df.detected_response == response) &
                                                 (credible_df.umi_count_tcr >= tcr_threshold)].copy() # & (credible_df.clonotype != 'None') & (credible_df.epitope != '0')
                            if exclude_clonotype_singlets:
                                sub_df.drop(sub_df[~sub_df.duplicated(subset='num_clonotype', keep=False)].index, inplace=True)
                            if exclude_specificity_singlets:
                                sub_df.drop(sub_df[~sub_df.duplicated(subset=['num_clonotype', 'epitope'], keep=False)].index, inplace=True)

                            gems = sub_df.gem.to_list()
                            mhc_read_count = sub_df.read_count_mhc.to_list()

                            xmax += len(np.unique(gems))
                            unique_gems.update(sub_df.gem.to_list())
                            unique_tcrs.update(sub_df.clonotype.to_list())

                            epitopes = sub_df.peptide_HLA.to_list()
                            colors = [epitope_to_color[ep] for ep in epitopes]

                            # How to show number of GEMs? Plot a different symbol if only one GEM?
                            #ax.scatter(gems, [-5]*len(gems))
                            scatter = ax.scatter(gems, mhc_read_count, marker=marker, c=colors) #, edgecolors='face', cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems, edgecolors='face', alpha=0.3

                        if i % 2 == 0:
                            plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)

                        xmin = xmax-0.5

                    from matplotlib.lines import Line2D

                    legend_elements = []
                    for epitope_label in all_epitopes:
                        legend_elements += [Line2D([0], [0], marker='o', color='w', label=epitope_label, markerfacecolor=epitope_to_color[epitope_label], markersize=10)]

                    legend1 = ax.legend(handles=legend_elements, ncol=7, loc=2, bbox_to_anchor=(0.02, -0.03))
                    ax.add_artist(legend1)

                    plt.tick_params(labelbottom=False, labelright=True, labelsize=8) #labelbottom=False, 
                    plt.xticks(rotation=90, size=2)

                    for marker, response in [('o', True), ('+', False)]:
                        plt.scatter([], [], c='k', label=str(response), marker=marker)
                    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Detected response', loc='upper left')

                    # Criteria
                    textstr = '\n'.join((
                        "Criteria",
                        "Min. barcode read count \t %i".expandtabs() %read_threshold,
                        "Min. TCR read count \t\t    %i".expandtabs() %tcr_threshold,
                        "Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
                        "Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
                    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
                    ax.text(0.12, 0.98, textstr,
                            fontsize=12,
                            horizontalalignment='left',
                            verticalalignment='top',
                            clip_on=False,
                            transform=ax.transAxes,
                            bbox=props)

                    plt.xlabel("%i GEMs (sectioned per clonotype (%i))" %(len(unique_gems), len(unique_tcrs)), fontsize=16)
                    plt.ylabel("pMHC barcode read counts", fontsize=16)
                    plt.title("MHC barcode read counts per GEM per clonotype", fontsize=20)
                    
                    if show:
                        plt.show()
                        print("OBS! Figures at not saved!")
                        return
                    if save_tuba:
                        plt.savefig(save_tuba + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    if save_sund:
                        plt.savefig(save_sund + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    plt.cla()   # Clear axis
                    plt.clf()   # Clear figure
                    plt.close(fig) # Close a figure window


# ## Peptide assayed

# In[10]:


def mhc_read_count_per_clonotype_peptide_assayed(credible_df, version=1, show=True, save_tuba=False, save_sund=False):
    import matplotlib as mpl
    mpl.rcParams['axes.grid.axis'] = 'y'#both
    
    # Epitope to color
    all_epitopes = credible_df.peptide_HLA.unique()
    col_epitopes = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'] * len(all_epitopes)
    epitope_to_color = dict()
    for i, epitope in enumerate(all_epitopes):
        epitope_to_color[epitope] = col_epitopes[i]
        
    # Detected response
    sortby = 'read_count_mhc' #'umi_count_tcr' #
    credible_df.sort_values(by=['num_clonotype', sortby], inplace=True)

    project = "read_count_per_clonotype_peptide_assayed/"
    
    for read_threshold in [1, 10]: # , 20, 50
        for tcr_threshold in [1, 10]: #, 20, 50
            for exclude_clonotype_singlets in [False, True]:
                for exclude_specificity_singlets in [False, True]:
                    #print(str(read_threshold), str(tcr_threshold), str(exclude_clonotype_singlets), str(exclude_specificity_singlets))

                    unique_gems, unique_tcrs = set(), set()

                    fig, ax = plt.subplots(figsize=(20, 10))

                    xmin, xmax = -0.5, 0

                    for i, clonotype in enumerate(credible_df.num_clonotype.unique()):
                        for marker, response in [('o', True), ('+', False)]:
                            sub_df = credible_df[(credible_df.num_clonotype == clonotype) &
                                                 (credible_df.read_count_mhc >= read_threshold) &
                                                 (credible_df.peptide_assayed == response) &
                                                 (credible_df.umi_count_tcr >= tcr_threshold)].copy() # & (credible_df.clonotype != 'None') & (credible_df.epitope != '0')
                            if exclude_clonotype_singlets:
                                sub_df.drop(sub_df[~sub_df.duplicated(subset='num_clonotype', keep=False)].index, inplace=True)
                            if exclude_specificity_singlets:
                                sub_df.drop(sub_df[~sub_df.duplicated(subset=['num_clonotype', 'epitope'], keep=False)].index, inplace=True)


                            gems = sub_df.gem.to_list()
                            mhc_read_count = sub_df.read_count_mhc.to_list()

                            xmax += len(np.unique(gems))
                            unique_gems.update(sub_df.gem.to_list())
                            unique_tcrs.update(sub_df.clonotype.to_list())

                            epitopes = sub_df.peptide_HLA.to_list()
                            colors = [epitope_to_color[ep] for ep in epitopes]
                            #tcr_umis = credible_df.umi_count_tcr.to_list()

                            # How to show number of GEMs? Plot a different symbol if only one GEM?
                            #ax.scatter(gems, [-5]*len(gems))
                            scatter = ax.scatter(gems, mhc_read_count, marker=marker, c=colors) #s=tcr_umis, , edgecolors='face', cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems, edgecolors='face', alpha=0.3

                        if i % 2 == 0:
                            plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)

                        xmin = xmax-0.5

                    from matplotlib.lines import Line2D

                    legend_elements = []
                    for epitope_label in all_epitopes:
                        legend_elements += [Line2D([0], [0], marker='o', color='w', label=epitope_label, markerfacecolor=epitope_to_color[epitope_label], markersize=10)]

                    legend1 = ax.legend(handles=legend_elements, ncol=7, loc=2, bbox_to_anchor=(0.02, -0.03))
                    ax.add_artist(legend1)

                    plt.tick_params(labelbottom=False, labelright=True, labelsize=8) #labelbottom=False, 
                    plt.xticks(rotation=90, size=2)

                    # OBS! DO NOT DELETE
                    #for size in [10, 50, 100]:
                    #    plt.scatter([], [], c='k', s=size, label=str(size) + " TCR UMIs", marker='o')
                    #legend2 = ax.legend(scatterpoints=1, frameon=False, labelspacing=1, title='TCR UMIs', loc='upper left')
                    #ax.add_artist(legend2)
                    #for marker, response in [('+', "no peptide assay")]: #('o', True), 
                    #    plt.scatter([], [], c='k', label=str(response), marker=marker)
                    #legend3 = ax.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Label indicator', loc='upper left')
                    #ax.add_artist(legend3)

                    for marker, response in [('o', True), ('+', False)]:
                        plt.scatter([], [], c='k', label=str(response), marker=marker)
                    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Peptide assay', loc='upper left')

                    # Criteria
                    textstr = '\n'.join((
                        "Criteria",
                        "Min. barcode read count \t %i".expandtabs() %read_threshold,
                        "Min. TCR read count \t\t    %i".expandtabs() %tcr_threshold,
                        "Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
                        "Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
                    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
                    ax.text(0.12, 0.98, textstr,
                            fontsize=12,
                            horizontalalignment='left',
                            verticalalignment='top',
                            clip_on=False,
                            transform=ax.transAxes,
                            bbox=props)

                    plt.xlabel("%i GEMs (sectioned per clonotype (%i))" %(len(unique_gems), len(unique_tcrs)), fontsize=16)
                    plt.ylabel("pMHC barcode read counts", fontsize=16)
                    plt.title("MHC barcode read counts per GEM per clonotype", fontsize=20)
                    
                    if show:
                        plt.show()
                        print("OBS! Figures at not saved!")
                        return
                    if save_tuba:
                        plt.savefig(save_tuba + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    if save_sund:
                        plt.savefig(save_sund + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    plt.cla()   # Clear axis
                    plt.clf()   # Clear figure
                    plt.close(fig) # Close a figure window


# # TCR read count per clonotype

# In[11]:


def tcr_read_count_per_clonotype_detected_response(credible_df, version=1, show=True, save_tuba=False, save_sund=False):
    mpl.rcParams['axes.grid.axis'] = 'y'#both
    
    # Epitope to color
    all_epitopes = credible_df.peptide_HLA.unique()
    col_epitopes = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'] * len(all_epitopes)
    epitope_to_color = dict()
    for i, epitope in enumerate(all_epitopes):
        epitope_to_color[epitope] = col_epitopes[i]
    
    sortby = 'umi_count_tcr'
    credible_df.sort_values(by=['num_clonotype', sortby], inplace=True)
    print(credible_df.shape)
    
    project = "tcr_read_count_per_clonotype_detected_response/"
    
    for read_threshold in [1, 2, 3]: #, 50
        for tcr_threshold in [1, 2, 3]: #, 20, 50
            for exclude_clonotype_singlets in [False, True]:
                for exclude_specificity_singlets in [False, True]:
                    unique_gems, unique_tcrs = set(), set()
                    
                    fig, ax = plt.subplots(figsize=(20, 10))

                    xmin, xmax = -0.5, 0

                    for i, clonotype in enumerate(credible_df.num_clonotype.unique()):
                        for marker, response in [('o', True)]: #, ('+', False)
                            sub_df = credible_df[(credible_df.num_clonotype == clonotype) &
                                                 (credible_df.umi_count_mhc >= read_threshold) &
                                                 (credible_df.umi_count_tcr >= tcr_threshold)].copy() #(credible_df.detected_response == response)
                            if exclude_clonotype_singlets:
                                sub_df.drop(sub_df[~sub_df.duplicated(subset='num_clonotype', keep=False)].index, inplace=True)
                            if exclude_specificity_singlets:
                                sub_df.drop(sub_df[~sub_df.duplicated(subset=['num_clonotype', 'epitope'], keep=False)].index, inplace=True)

                            gems = sub_df.gem.to_list()
                            tcr_read_count = sub_df.umi_count_tcr.values

                            xmax += len(np.unique(gems))
                            
                            unique_gems.update(sub_df.gem.to_list())
                            unique_tcrs.update(sub_df.clonotype.to_list())

                            epitopes = sub_df.peptide_HLA.to_list()
                            colors = [epitope_to_color[ep] for ep in epitopes]

                            # How to show number of GEMs? Plot a different symbol if only one GEM?
                            #ax.scatter(gems, [-5]*len(gems))
                            scatter = ax.scatter(gems, tcr_read_count, marker=marker, c=colors) #, edgecolors='face', cmap='viridis_r', norm=plt.Normalize(vmin=0, vmax=1), s=gems, edgecolors='face', alpha=0.3

                        if i % 2 == 0:
                            plt.axvspan(xmin, xmax-0.5, facecolor='0.7', alpha=0.1)

                        xmin = xmax-0.5

                    from matplotlib.lines import Line2D

                    legend_elements = []
                    for epitope_label in all_epitopes:
                        legend_elements += [Line2D([0], [0], marker='o', color='w', label=epitope_label, markerfacecolor=epitope_to_color[epitope_label], markersize=10)]

                    legend1 = ax.legend(handles=legend_elements, ncol=8, loc=2, bbox_to_anchor=(0.02, -0.03))
                    ax.add_artist(legend1)

                    plt.tick_params(labelbottom=False, labelright=True, labelsize=8) #labelbottom=False, 
                    plt.xticks(rotation=90, size=2)

                    for marker, response in [('o', True), ('+', False)]:
                        plt.scatter([], [], c='k', label=str(response), marker=marker)
                    plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Detected response', loc='upper right')

                    # Criteria
                    textstr = '\n'.join((
                        "Criteria",
                        "Min. BC UMI count \t %i".expandtabs() %read_threshold,
                        "Min. TCR UMI count \t\t    %i".expandtabs() %tcr_threshold,
                        "Exclude clonotype singlets \t %s".expandtabs() %str(exclude_clonotype_singlets),
                        "Exclude specificity singlets \t   %s".expandtabs() %str(exclude_specificity_singlets)))
                    props = dict(boxstyle='square', fc='white', ec='grey', alpha=0.5)
                    ax.text(0.05, 0.95, textstr,
                            fontsize=12,
                            horizontalalignment='left',
                            verticalalignment='top',
                            clip_on=False,
                            transform=ax.transAxes,
                            bbox=props)

                    plt.xlabel("%i GEMs (sectioned per clonotype (%i))" %(len(unique_gems), len(unique_tcrs)), fontsize=16)
                    plt.ylabel("TCR UMI counts", fontsize=16)
                    plt.title("TCR UMI counts per GEM per clonotype", fontsize=20)
                    
                    if show:
                        plt.show()
                        print("OBS! Figures at not saved!")
                        return
                    if save_tuba:
                    	# OBS! save_tuba + project + 
                        plt.savefig("v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(4, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    if save_sund:
                        plt.savefig(save_sund + project + "v%i.b%i.t%i.ecs_%s.ess_%s.pdf" %(version, read_threshold, tcr_threshold, exclude_clonotype_singlets, exclude_specificity_singlets), bbox_inches='tight')
                    plt.cla()   # Clear axis
                    plt.clf()   # Clear figure
                    plt.close(fig) # Close a figure window


# In[ ]:




