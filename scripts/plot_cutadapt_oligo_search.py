#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')


# # Args

# In[2]:


#EXP = "exp5"
#PLATFORM = "IONTORRENT"


# In[4]:


FIG_DIR = snakemake.params[0]


# # Distribution of reads annotated with building blocks

# In[5]:


#FILE = "/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC_" + PLATFORM + "/reports/cutadapt_oligo_search/headers/oligo_annotated_headers.uniq_count.lst"


# In[6]:


head_df = pd.read_csv(snakemake.input[0], sep=" ", names=("counts", "tso", "b_primer", "anneal", "cd8_primer", "mhc_primer"))


# In[7]:


#head_df.head()


# In[8]:


head_df['a_primer'] = np.where(head_df.cd8_primer == 'no_adapter', head_df.mhc_primer, head_df.cd8_primer)
head_df.replace('no_adapter', np.NaN, inplace=True)
head_df.a_primer = head_df.a_primer.str.slice_replace(start=7, repl='A')
head_df['match'] = head_df.loc[:, ['tso', 'b_primer', 'anneal', 'a_primer']].count(axis=1)


# In[9]:


head_df.replace(np.NaN, '', inplace=True)
head_df['name'] = head_df.tso + " " + head_df.b_primer + " " + head_df.anneal + " " + head_df.a_primer
head_df.loc[head_df.match == 0, 'name'] = "no adapter"


# In[10]:


#head_df


# In[11]:


title = "Number of construct annotations in reads (3,035,311 reads)"
total_sizes = list()
for i in range(5):
    total_sizes.append(head_df[head_df.match == i].counts.sum())
    
labels = tuple(zip(list(range(5)), total_sizes))

fig1, ax1 = plt.subplots()
wedges, text = ax1.pie(total_sizes, shadow=True, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.title(title)

ax1.legend(wedges, labels,
          title="Building blocks",
          loc="center left",
          bbox_to_anchor=(0.8, 0, 0.5, 1))

plt.savefig(FIG_DIR + "pie_all_present.pdf", bbox_inches='tight') 
#plt.show()

print(total_sizes)


# In[12]:



for i in range(1,4):
    title = "Distribution of constructs (when %i present)" %i
    sizes = head_df[head_df.match == i].groupby('name')['counts'].sum().to_list()
    labels = [" ".join(j.strip().split()) for j in head_df[head_df.match == i].groupby('name')['counts'].sum().index]
    fig1, ax1 = plt.subplots()
    wedges, text = ax1.pie(sizes, shadow=False, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title(title)

    ax1.legend(wedges, labels,
              title="Building blocks",
              loc="center left",
              bbox_to_anchor=(0.8, 0, 0.5, 1))

    plt.savefig(FIG_DIR + "pie_%i_present.pdf" %i, bbox_inches='tight')
    #plt.show()


# In[ ]:




