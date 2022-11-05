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

def get_unique_entries(df):
    return np.where((df.cdr3_TRA==a) & (df.cdr3_TRB==b), False, True)

def get_sample_size():
    inter_entries = get_unique_entries(inter_chains)
    inter_indexes = inter_chains[inter_entries].index.to_list()
    return min(sum(get_unique_entries(group)), len(inter_indexes))

def add_similarity_scores(ai,bi):
    mat_a = sim_tra.loc[ai, a].reset_index(drop=True).T.reset_index(drop=True).T
    mat_b = sim_trb.loc[bi, b].reset_index(drop=True).T.reset_index(drop=True).T 
    return mat_a.add(mat_b)

def get_intra_similarity(cdr3_TRAs, cdr3_TRBs):
    unique_entries = get_unique_entries(group)
    unique_entry_indexes = group[unique_entries].index.to_list()

    sample_size = get_sample_size()

    intra_sample = random.sample(unique_entry_indexes, sample_size)
    intra_a = group.loc[intra_sample, 'cdr3_TRA'].values
    intra_b = group.loc[intra_sample, 'cdr3_TRB'].values
    combined_similarity = add_similarity_scores(intra_a, intra_b)
    return {'score': combined_similarity.max(),
            'fraction': sum(combined_similarity > 1.8)/len(combined_similarity)}

def get_inter_similarity(cdr3_TRAs, cdr3_TRBs):
    # OBS! make sure the size to sample from matches the number og unique entries intra_similarity! 
    unique_entries = get_unique_entries(inter_chains)
    unique_entry_indexes = inter_chains[unique_entries].index.to_list()

    sample_size = get_sample_size()

    inter_sample = random.sample(unique_entry_indexes, sample_size)
    inter_a = inter_chains.loc[inter_sample, 'cdr3_TRA'].values
    inter_b = inter_chains.loc[inter_sample, 'cdr3_TRB'].values
    combined_similarity = add_similarity_scores(inter_a, inter_b) 
    return {'score': combined_similarity.max(),
            'fraction': sum(combined_similarity > 1.8)/len(combined_similarity)}

def paired_t_test(x1 ,x2):
    assert len(x1) == len(x2)
    statistic, pvalue = stats.ttest_rel(x1, x2)
    if (pvalue/2.0 < 0.05) & (statistic > 0) & (len(x1) > 9):
        return {'test':True, 'pvalue':pvalue}
    else:
        return {'test':False, 'pvalue':pvalue}

def add_number_of_observations(intra_lst, inter_lst):
    for box, lst in enumerate([intra_lst, inter_lst], start=1):
        if lst:
            median = statistics.median(lst)
            plt.text(box, median, "n: %i" %len(lst), ha='center', va='bottom')

def add_significance_bar(intra_lst, inter_lst):
    t = paired_t_test(intra_lst, inter_lst)#['pvalue']
    if t['test'] and t['pvalue'] < 0.05:
        pass
    else:
        return
    
    y0 =  max(max(intra_lst), max(inter_lst))
    y1 = y0 * 1.02
    y2 = y0 * 1.025
    y3 = y0 * 1.03
    y4 = y0 * 1.035

    plt.plot([1,1,2,2], [y1,y2,y2,y1], lw=1.5, c='k')
    plt.text(1.5, y3, "p = %.6f" %t['pvalue'], ha='center', va='bottom', color='k')
    plt.plot(1, y4)

def plot_boxplot(intra_lst, inter_lst):
    plt.figure(figsize=(3,5))
    plt.boxplot([intra_lst, inter_lst], labels=['intra', 'inter'], widths=(0.5, 0.5))
    plt.title('Pooled')
    plt.xlim(0.6, 2.4)
    plt.ylabel("Similarity")
    
    add_number_of_observations(intra_lst, inter_lst)
    add_significance_bar(intra_lst, inter_lst)

    plt.savefig(get_file('boxplot'), bbox_inches='tight')
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close()

def plot_pieplot(significant_count, total_peptides):
    plt.pie([significant_count, total_peptides-significant_count],
            labels=['significant', 'insignificant'],
            autopct=lambda p: '{:.0f} ({:.0f}%)'.format(p * total_peptides / 100, p))
    plt.title("Proportion of significant outcomes (%i)" %total_peptides)
    plt.savefig(get_file('pieplot'), bbox_inches='tight')
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close()


# ARGS
min_binding_concordance = 0.05
n_simulations = 100

#if __name__ == '__main__':
## ## Input
INPUT_DF = snakemake.input.df
INPUT_GEM = snakemake.input.unique_gems
INPUT_PEP = snakemake.input.peptides

INPUT_TRA = snakemake.input.sim_tra
INPUT_TRB = snakemake.input.sim_trb

OUTPUT_X = snakemake.output.file_X
OUTPUT_Y = snakemake.output.file_Y
OUTPUT_Zs = snakemake.output.siml_Z
OUTPUT_Zg = snakemake.output.gems_Z

# # Load
df_input = pd.read_csv(INPUT_DF, low_memory=False)
unique_gem_file = iter(sort_input_files(list(set(INPUT_GEM))))

input_peptides = pd.read_excel(INPUT_PEP)
total_peptides = len(input_peptides.Sequence.unique()) - 1 # Perhaps choose a selected of expected peptides (peptides we expect to be recognized by TCR) --> ranked peptides in annotation file?

sim_tra = pd.read_csv(INPUT_TRA, index_col=0).rename(index={'missing':''}, columns={'missing':''})
sim_trb = pd.read_csv(INPUT_TRB, index_col=0).rename(index={'missing':''}, columns={'missing':''})

sim_tra = sim_tra[~sim_tra.index.duplicated()]
sim_trb = sim_trb[~sim_trb.index.duplicated()]

# # Main
bc_umis, tcr_umis = extract_minimum_umi_counts(INPUT_GEM)

# Initialize comparison arrays
X, Y = np.meshgrid(bc_umis, tcr_umis) #sorted(bc_umis, reverse=True), sorted(tcr_umis, reverse=True)
z_similarity_matrix = np.zeros((len(tcr_umis), len(bc_umis)))
z_gem_matrix = np.zeros((len(tcr_umis), len(bc_umis)))

# Setup for comparison
for x, min_mhc_umi_count in enumerate(bc_umis):
    for y, min_tcr_umi_count in enumerate(tcr_umis):
        print('b: {}, t: {}'.format(min_mhc_umi_count, min_tcr_umi_count))
        filename = next(unique_gem_file)
        unique_gems = np.loadtxt(filename, dtype='U20')
        print(filename)

        # Condition dataframe
        df = df_input[df_input.gem.isin(unique_gems)]
        df = df.replace('unknown', np.nan).dropna(subset=['cdr3_TRA', 'cdr3_TRB'])

        intra_score = list()
        inter_score = list()

        significant_simulations = list()
        for simulation in range(n_simulations):
            significant_count = 0
            for peptide, group in df.groupby('peptide_HLA'):
                if len(group) == 1:
                    continue
                if len(group.drop_duplicates(['cdr3_TRA','cdr3_TRB'])) == 1:
                    continue
                    
                inter_chains = df.loc[df.peptide_HLA != peptide, ['cdr3_TRA', 'cdr3_TRB']]

                intra_score_peptide = list()
                inter_score_peptide = list()  
                
                cdr3_TRAs = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRA.values #.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB'])
                cdr3_TRBs = group.drop_duplicates(subset=['cdr3_TRA','cdr3_TRB']).cdr3_TRB.values
                
                for index, (a,b) in enumerate(zip(cdr3_TRAs, cdr3_TRBs)):  
                    intra = get_intra_similarity(cdr3_TRAs, cdr3_TRBs)
                    inter = get_inter_similarity(cdr3_TRAs, cdr3_TRBs)
                    
                    intra_score_peptide.append(intra['score'])
                    inter_score_peptide.append(inter['score'])
                    
                    intra_score.append(intra['score'])
                    inter_score.append(inter['score'])
                    
                # Statistics
                # Calc p-value for each peptide-plateau (count number of intra significant peptide-plateaus)
                # Count fraction of GEMs above 1.8 similarity
                if paired_t_test(intra_score_peptide, inter_score_peptide)['test']:
                    significant_count += 1
            significant_simulations.append(significant_count)
        significant_count = statistics.median(significant_simulations)
        print(significant_count)
        plot_pieplot(significant_count, total_peptides)

        if (min_mhc_umi_count == max(bc_umis)) & (min_tcr_umi_count == max(tcr_umis)):
            tnobs = len(intra_score)
        else:
            intra_score = random.sample(intra_score, tnobs)
            inter_score = random.sample(inter_score, tnobs)
        plot_boxplot(intra_score, inter_score)

        z_similarity_matrix[y,x] = significant_count/float(total_peptides)
        z_gem_matrix[y,x] = len(unique_gems)

np.savetxt(OUTPUT_X , X, fmt='%d')
np.savetxt(OUTPUT_Y , Y, fmt='%d')
np.savetxt(OUTPUT_Zs, z_similarity_matrix)
np.savetxt(OUTPUT_Zg, z_gem_matrix)

