#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')

fig, ax = plt.subplots(len(snakemake.input),1, sharex=True, figsize=(10,10))
#plt.subplots_adjust(hspace=0.5)
for i, filename in enumerate(snakemake.input):
	df = pd.read_csv(snakemake.input[i], r"\s+", names=['frequency', 'read_length'])

	# Find median
	df.sort_values('read_length', inplace=True)
	cumsum = df.frequency.cumsum()
	cutoff = df.frequency.sum() / 2.0
	median = df.read_length[cumsum >= cutoff].iloc[0]

	fn = filename.split("/")[-1].split(".")[0]
	if fn == 'cut_adapter':
		title = 'After adapter trimming'
	elif fn == 'longranger':
		title = 'After GEM annotation'
	elif fn == 'longranger_clean':
		title = 'After filtering reads without GEM annotation'
	elif fn == 'output':
		title = 'After mapping to BC library'
	else:
		title = fn

	ax[i].bar(df.read_length.values, df.frequency.values)
	ax[i].set_title(title, fontsize=10)
	ax[i].axvline(median, color='k', linestyle='dashed', linewidth=1)

fig.text(0.5, 0.04, 'Read lengths', ha='center')
fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
fig.suptitle('Read length distributions')

plt.savefig(snakemake.output[0], bbox_inches='tight')
plt.show()
