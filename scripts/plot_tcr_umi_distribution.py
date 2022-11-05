#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')

def get_bins(x):
	if len(x) == 0:
		print(len(x))
		return 1
	q25, q75 = np.percentile(x,[.25,.75])
	if q25 == q75:
		print(x.max()+2)
		return int(x.max()+2)
	bin_width = 2*(q75 - q25)*len(x)**(-1/3)
	return round((x.max() - x.min())/bin_width)

fig, ax = plt.subplots(len(snakemake.input),1, sharex=True, figsize=(10,10))
#plt.subplots_adjust(hspace=0.5)
for i, filename in enumerate(snakemake.input):
	fn = filename.split("/")[-1].rsplit(".", 1)[0]
	print(fn)
	if fn == 'umi_per_gem.raw':
		df = pd.read_csv(snakemake.input[i], sep=r'\s+', header=None, names=['gem', 'freq'])
		title = 'Raw TCR reads'
		x = df.freq
		a,b = [],[]
		upper_limit = max(x)
	elif fn == 'umi_per_gem.filtered':
		df = pd.read_csv(snakemake.input[i], sep=r'\s+', header=None, names=['gem', 'freq'])
		title = 'TCR reads filtered away by cellranger'
		x = df.freq
		a,b = [],[]

	df = pd.read_csv(snakemake.input[i])
	if fn == 'all_contig_annotations':
		title = 'All TCR contigs'
		a = df[df.chain == 'TRA'].umis.dropna()
		b = df[df.chain == 'TRB'].umis.dropna()
		x = []
	elif fn == 'tcr.clean.augmented':
		title = 'Selected TCR contigs'
		a = df.umi_count_TRA.dropna()
		b = df.umi_count_TRB.dropna()
		x = []
	elif fn == 'tcr_barcode.cleaned':
		title = 'TCR contigs also containing pMHC barcodes'
		a = df.umi_count_TRA.dropna()
		b = df.umi_count_TRB.dropna()
		x = []
	#else:
	#	title = fn
	#	print(title)

	ax[i].hist(x, bins=get_bins(x), alpha=0.7, label='Any')
	ax[i].hist(b, bins=get_bins(b), alpha=0.7, label='TRB')
	ax[i].hist(a, bins=get_bins(a), alpha=0.7, label='TRA')
	ax[i].legend()
	ax[i].set_title(title, fontsize=10)
	ax[i].set_xscale('log')
	ax[i].set_xlim(0.91, upper_limit)

fig.text(0.5, 0.04, 'UMIs per GEM', ha='center')
fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical')
fig.suptitle('UMI distributions for GEMs containing TCRs')
#fig.set_xscale('log')
#fig.set_xlim(0.9, upper_limit)

plt.savefig(snakemake.output[0], bbox_inches='tight')
plt.show()
