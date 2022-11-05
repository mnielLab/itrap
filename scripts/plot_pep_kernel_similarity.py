#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import random


# In[2]:


plt.style.use('ggplot')


# # Args

# In[4]:


EXP = "exp3"
PLATFORM = "IONTORRENT"


# In[5]:


MAPPING = 'KMA' # BLAST
BARCODE_SYSTEM = 'AKB' #'AKB' #


# In[6]:


def pd_fill_diagonal(df_matrix, value=0): 
    mat = df_matrix.values
    n = mat.shape[0]
    mat[range(n), range(n)] = value
    return pd.DataFrame(mat)


# In[7]:


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
    if pvalue < 0.05:
        y, h = simila.max().max() + simila.max().max() * 0.01, 2 * simila.max().max() * 0.005
        plt.plot([1, 1, 2, 2], [y, y+h, y+h, y], lw=1.5, c='k')
        plt.text(1.5, y+h, "p = %.6f" %pvalue, ha='center', va='bottom', color='k')
    plt.savefig(OUT_DIR + "overall/min%i/%s.pdf" %(min_read_count, title), bbox_inches='tight')
    plt.show()


# ## Input

# In[8]:


IN_FILE = ("/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" +
           EXP + "_CAT_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
           "/similarity_assessment/peptide.csv")


# ## Output dir

# In[9]:


OUT_DIR = ("/Volumes/tuba/herpov/tcr-pmhc-sc-project/plots/" +
           EXP + "_" + PLATFORM + "_" + MAPPING + "_" + BARCODE_SYSTEM +
           "/similarity_assessment/peptide/")


# # Main

# In[10]:


total_df = pd.read_csv(IN_FILE)


# In[12]:


total_df.rename(columns={'peptide_x':'peptide'}, inplace=True)


# ## Plot heatmap

# In[14]:


color_cdr3 = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'] * total_df.clonotype.unique().shape[0]
my_palette = dict(zip(total_df.clonotype.unique(), color_cdr3))
row_colors = total_df.clonotype.map(my_palette)


# In[13]:


# https://python-graph-gallery.com/404-dendrogram-with-heat-map/
# https://stackoverflow.com/questions/27988846/how-to-express-classes-on-the-axis-of-a-heatmap-in-seaborn
sns.clustermap(total_df.iloc[:, 7:], standard_scale=1,
               row_cluster=False, col_cluster=False,
               row_colors=row_colors,
               xticklabels=False, yticklabels=False, cmap='seismic')# col_colors=row_colors
plt.savefig(OUT_DIR + "heatmap/all.png")
plt.show()


# # V3

# In[16]:


mhc_read_counts = [1, 10]
min_binding_concordance = 0.4
for min_read_count in sorted(mhc_read_counts, reverse=True):
    intra_scores = list()
    inter_scores = list()
    significance_count = 0

    df = total_df[total_df.read_counts_mhc >= min_read_count].set_index('gem') #NB! Should I remove GEMs with same peptides and clonotypes?

    clonotypes = total_df.clonotype.unique()
    for clonotype in clonotypes:
        # Intra
        i = df[(df.clonotype == clonotype) & (df.binding_concordance >= min_binding_concordance)].drop_duplicates(subset=['peptide']).index
        matrix = pd_fill_diagonal(df.loc[i,i], np.nan) #.iloc[i,i+5]
        matrix.columns = i
        intra = matrix.max()

        # Inter
        indexes = df[(~df.clonotype.isin([clonotype])) & (df.binding_concordance >= min_binding_concordance)].drop_duplicates(subset=['peptide']).index
        j = random.sample(indexes.to_list(), len(i))
        matrix = df.loc[j, i]
        inter = matrix.max()

        simila = pd.DataFrame([intra, inter])
        labels = ['intra', 'inter']

        nobs = simila.shape[1]

        if nobs > 1:
            # Make boxplot
            fig = plt.figure(figsize=(3,5))
            plt.boxplot(simila, widths=(0.5, 0.5))
            plt.xticks(range(1, len(labels) + 1), labels)
            plt.ylabel("Similarity")
            plt.title("%s" %clonotype) #Inter vs. intra peptide similarity for 
            plt.xlim(0.6,2.4)

            # Text on boxplot
            medians = simila.median(axis=1).values
            for l in range(1, len(labels)+1):
                plt.text(l, medians[l-1], "n: %i" %nobs, ha='center', va='bottom')#+medians[l-1]*0.005

            # Statistical annotation
            statistic, pvalue = stats.ttest_ind(simila.iloc[[0],:].values[0],
                                                simila.iloc[[1],:].values[0],
                                                equal_var=False)
            if pvalue < 0.05:
                significance_count += 1

                y, h = simila.max().max() + simila.max().max() * 0.01, 2 * simila.max().max() * 0.005
                plt.plot([1, 1, 2, 2], [y, y+h, y+h, y], lw=1.5, c='k')
                plt.text(1.5, y+h, "*", ha='center', va='bottom', color='k')    

            plt.savefig(OUT_DIR + "individual/min%i/%s.pdf" %(min_read_count, clonotype), bbox_inches='tight')
            #plt.show()
            plt.cla()   # Clear axis
            plt.clf()   # Clear figure
            plt.close(fig)
            
            # Store similarity scores for peptides with minimum 9 TCRs
            intra_scores += intra.to_list()
            inter_scores += inter.to_list()

    plt.pie([significance_count, len(clonotypes)-significance_count],
            labels=['significant difference', 'insignificant'],
            autopct=lambda p: '{:.0f}'.format(p * len(clonotypes) / 100))
    plt.title("Distribution of significant outcomes (%i)" %len(clonotypes))
    plt.savefig(OUT_DIR + "overall/min%i/significant_outcomes.pdf" %min_read_count, bbox_inches='tight')
    plt.show()
    
    if min_read_count == max(mhc_read_counts):
        tnobs = len(intra_scores)
    else:
        intra_scores = random.sample(intra_scores, tnobs)
        inter_scores = random.sample(inter_scores, tnobs)
    
    plot_intra_vs_inter(pd.DataFrame({'intra': intra_scores, 'inter': inter_scores}, columns=['intra', 'inter']).T,
                        "all")

