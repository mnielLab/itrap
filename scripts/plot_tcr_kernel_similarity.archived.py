#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import random
import re
import statistics

plt.style.use('ggplot')

## # Args
#EXP = "exp6"
#PLATFORM = "IONTORRENT"
#MAPPING = 'KMA' # BLAST
#BARCODE_SYSTEM = 'AKB' #'AKB' #

def get_file(sub_level):
    #print('\/%s\/b%i\.t%i\.'%(sub_level, min_mhc_umi_count, min_tcr_umi_count))
    for filename in snakemake.output.plot:
        if re.search('\/%s\/b%i\.t%i\.'%(sub_level, min_mhc_umi_count, min_tcr_umi_count), filename):
            #print(filename)
            return filename

def sort_input_files(input_files):
	def atoi(text):
	    return int(text) if text.isdigit() else text

	def natural_keys(text):
	    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

	input_files.sort(key=natural_keys, reverse=True)

	return input_files

def extract_minimum_umi_counts(input_files):
	def do_search(obj):
		return re.search('b(\d+).t(\d+)', obj)

	def get_int(obj):
		return int(do_search(obj).group(1)), int(do_search(obj).group(2))

	def list_min_umi(input_files):
		umi_count = [get_int(filename) for filename in input_files]
		return list(zip(*umi_count))

	def rev_sort_umi(input_files):
		return [sorted(list(set(l)), reverse=True) for l in list_min_umi(input_files)]

	return rev_sort_umi(input_files)[0], rev_sort_umi(input_files)[1]

def reduce_df_by_unique_gems(total_df, filename):
    #for filename in snakemake.input.unique_gems:
    if re.search('b%i.t%i.'%(min_mhc_umi_count, min_tcr_umi_count), filename):
    	unique_gems = np.loadtxt(filename, dtype='U20')
    	total_df = total_df.loc[total_df.gem.isin(unique_gems), np.append(total_df.columns[:12], total_df.columns.intersection(unique_gems))]
    	total_df = total_df.set_index('gem')

    	assert total_df.equals(total_df[(total_df.umi_count_mhc >= min_mhc_umi_count) & (total_df.umi_count_tcr >= min_tcr_umi_count)])
    	return total_df

def pd_fill_diagonal(df_matrix, value=0): 
    mat = df_matrix.values
    n = mat.shape[0]
    mat[range(n), range(n)] = value
    return pd.DataFrame(mat)

def plot_significance_bar():
	b, h = simila.max().max() + simila.max().max() * 0.01, 2 * simila.max().max() * 0.005
	plt.plot([1, 1, 2, 2], [b, b+h, b+h, b], lw=1.5, c='k')
	plt.text(1.5, b+h, "*", ha='center', va='bottom', color='k')

def plot_individual_boxplots(): # OBS! get_file() is changed!!
    # Make boxplot
    plt.figure(figsize=(3,5))
    plt.boxplot(simila, widths=(0.5, 0.5))
    plt.xticks(range(1, len(labels) + 1), labels)
    plt.ylabel("Similarity")
    plt.title("%s" %pep) #Inter vs. intra peptide similarity for 
    plt.xlim(0.6,2.4)

    # Text on boxplot
    medians = simila.median(axis=1).values
    for l in range(1, len(labels)+1):
        plt.text(l, medians[l-1], "n: %i" %nobs, ha='center', va='bottom')#+medians[l-1]*0.005

    # Statistical annotation OBS! Statistics are calculated outside the function in 'main'
    #statistic, pvalue = stats.ttest_ind(simila.iloc[[0],:].values[0],
    #                                    simila.iloc[[1],:].values[0],
    #                                    equal_var=False)
    if (pvalue/2.0 < 0.05) & (statistic > 0): # One-tailed 'greater-than' test
        plot_significance_bar()    
    #plt.savefig(OUT_DIR + "individual/min%i/%s.pdf" %(min_mhc_umi_count, pep), bbox_inches='tight')
    plt.savefig('%s%s.pdf'%(get_file('individual'), pep), bbox_inches='tight')
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close()

def plot_significant_outcomes():
    plt.pie([significance_count, len(peptides)-significance_count],
            labels=['significant difference', 'insignificant'],
            autopct=lambda p: '{:.0f}'.format(p * len(peptides) / 100))
    plt.title("Distribution of significant outcomes (%i)" %len(peptides))
    #plt.savefig(OUT_DIR + "pooled/min%i/significant_outcomes.pdf" %min_mhc_umi_count, bbox_inches='tight')
    #plt.savefig('%s%s.pdf'%(get_file('pooled'), 'significant_outcomes'), bbox_inches='tight')
    plt.savefig(get_file('significant_outcomes'), bbox_inches='tight')
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close()

def plot_intra_vs_inter(simila, title):
    labels = ['intra', 'inter']

    plt.figure(figsize=(3,5))
    plt.boxplot(simila, widths=(0.5, 0.5))
    plt.xticks(range(1, len(labels) + 1), labels)
    plt.ylabel("Similarity")
    plt.title("%s" %title)
    plt.xlim(0.6,2.4)

    # Text
    medians = simila.median(axis=1).values
    nobs = simila.shape[1]
    for l in range(1, len(labels)+1):
        plt.text(l, medians[l-1], "n: %i" %nobs, ha='center', va='bottom')#+medians[l-1]*0.005

    # Statistical annotation
    statistic, pvalue = stats.ttest_ind(simila.iloc[[0],:].values[0],
                                        simila.iloc[[1],:].values[0],
                                        equal_var=False)
    if (pvalue/2.0 < 0.05) & (statistic > 0): # One-tailed 'greater-than' test
        y, h = simila.max().max() + simila.max().max() * 0.01, 2 * simila.max().max() * 0.005
        plt.plot([1, 1, 2, 2], [y, y+h, y+h, y], lw=1.5, c='k')
        plt.text(1.5, y+h, "p = %.6f" %pvalue, ha='center', va='bottom', color='k')
    #plt.savefig(OUT_DIR + "pooled/min%i/%s.pdf" %(min_mhc_umi_count, boxplot), bbox_inches='tight')
    plt.savefig(get_file('boxplot'), bbox_inches='tight')
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close()



## ## Input
#IN_FILE = ("/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" +
#           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
#           "/similarity_assessment/cdr3.csv")
#
#GEM_FIL = ("/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" +
#           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
#           "/specificity_matrix/peptide_per_clonotype_by_gem_size/num_clonotype/exclude_single-chain_TCRs/unique_gems/b2.t2.ecs_False.ess_False.lst")
#
## ## Output dir
#OUT_DIR = ("/Volumes/tuba/herpov/tcr-pmhc-sc-project/plots/" +
#           EXP + "_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
#           "/similarity_assessment/cdr3/")

# # Main
total_df = pd.read_csv(snakemake.input.df, low_memory=False)
peptides = total_df.peptide.unique() # Now peptides contain all barcode labeled peptides
unique_gem_file = iter(sort_input_files(list(set(snakemake.input.unique_gems))))


#gem_files = set(snakemake.input.unique_gems)
#for GEM_FIL in snakemake.input.unique_gems:
#    unique_gems = np.loadtxt(GEM_FIL, dtype='U20')
#    total_df = total_df.loc[total_df.gem.isin(unique_gems), np.append(total_df.columns[:12], total_df.columns.intersection(unique_gems))] # OBS! What number of columns?!
    
    ## ## Plot heatmap
    #color_cdr3 = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'] * total_df.peptide.unique().shape[0]
    #my_palette = dict(zip(total_df.peptide.unique(), color_cdr3))
    #row_colors = total_df.peptide.map(my_palette)
    ## https://python-graph-gallery.com/404-dendrogram-with-heat-map/
    ## https://stackoverflow.com/questions/27988846/how-to-express-classes-on-the-axis-of-a-heatmap-in-seaborn
    #sns.clustermap(total_df.iloc[:, 10:], standard_scale=1,
    #               row_cluster=False, col_cluster=False,
    #               row_colors=row_colors,
    #               xticklabels=False, yticklabels=False, cmap='seismic')# col_colors=row_colors
    ##plt.savefig(OUT_DIR + "heatmap/all.png")

bc_umis, tcr_umis = extract_minimum_umi_counts(snakemake.input.unique_gems)
# Initialize comparison arrays
#bc_umis = [1, 2, 3, 4, 5]
#tcr_umis = [1, 2, 4, 6, 8, 10]
X, Y = np.meshgrid(bc_umis, tcr_umis) #sorted(bc_umis, reverse=True), sorted(tcr_umis, reverse=True)
z_similarity_matrix = np.zeros((len(tcr_umis), len(bc_umis)))
z_gem_matrix = np.zeros((len(tcr_umis), len(bc_umis)))

min_binding_concordance = 0.05
for x, min_mhc_umi_count in enumerate(bc_umis): #sorted(bc_umis, reverse=True) #sorted(tcr_umis, reverse=True)
    for y, min_tcr_umi_count in enumerate(tcr_umis):
        filename = next(unique_gem_file)
        unique_gems = np.loadtxt(filename, dtype='U20')
        df = total_df.loc[total_df.gem.isin(unique_gems), np.append(total_df.columns[:12], total_df.columns.intersection(unique_gems))]
        df = df.set_index('gem')
        print(df.columns.to_list())
        assert df.equals(df[(df.umi_count_mhc >= min_mhc_umi_count) & (df.umi_count_tcr >= min_tcr_umi_count)]), print(df.umi_count_mhc.min(), min_mhc_umi_count, df.umi_count_tcr.min(), min_tcr_umi_count)

        intra_scores = list()
        inter_scores = list()

        n_simulations = 100
        significance_simulations = list()
        for simulation in range(n_simulations):
            significance_count = 0
            for pep in peptides:
    			# Intra
                i = df[(df.peptide == pep) & (df.binding_concordance >= min_binding_concordance)].drop_duplicates(subset=['cdr3']).index
                matrix = pd_fill_diagonal(df.loc[i,i], np.nan) #.iloc[i,i+5]
                matrix.columns = i
                intra = matrix.max()

                # Inter
                indexes = df[(~df.peptide.isin([pep, 'SLAAYIPRL'])) & (df.binding_concordance >= min_binding_concordance)].drop_duplicates(subset=['cdr3']).index
                j = random.sample(indexes.to_list(), len(i)) # should it be len(i)-1 ?!
                matrix = df.loc[j, i]
                inter = matrix.max()

                simila = pd.DataFrame([intra, inter])
                labels = ['intra', 'inter']

                nobs = simila.shape[1]
                if nobs > 1:
                    # Should be a paired t-test?
                    statistic, pvalue = stats.ttest_ind(simila.iloc[[0],:].values[0], simila.iloc[[1],:].values[0], equal_var=False)
                    # One-tailed 'greater-than' test
                    if (pvalue/2.0 < 0.05) & (statistic > 0):
                        significance_count += 1

                    #plot_individual_boxplots()
                    # Store similarity scores for peptides with minimum 9 TCRs
                    intra_scores += intra.to_list()
                    inter_scores += inter.to_list()

            significance_simulations.append(significance_count)
        assert len(significance_simulations) == n_simulations
        print(significance_simulations)
        significance_count = statistics.median(significance_simulations)

        plot_significant_outcomes()

        z_similarity_matrix[y,x] = significance_count/float(len(peptides))
        z_gem_matrix[y,x] = len(unique_gems)
        assert len(unique_gems) == df.shape[0]

        if (min_mhc_umi_count == max(bc_umis)) & (min_tcr_umi_count == max(tcr_umis)):
            tnobs = len(intra_scores)
        else:
            intra_scores = random.sample(intra_scores, tnobs)
            inter_scores = random.sample(inter_scores, tnobs)

        plot_intra_vs_inter(pd.DataFrame({'intra': intra_scores, 'inter': inter_scores}, columns=['intra', 'inter']).T, "Pooled")

np.savetxt(snakemake.output.file_X, X, fmt='%d')
np.savetxt(snakemake.output.file_Y, Y, fmt='%d')
np.savetxt(snakemake.output.siml_Z, z_similarity_matrix)
np.savetxt(snakemake.output.gems_Z, z_gem_matrix)

