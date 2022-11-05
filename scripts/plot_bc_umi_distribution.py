#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')

def get_bins(x):
	if len(x) == 0:
		return 1
	q25, q75 = np.percentile(x,[.25,.75])
	if q25 == q75:
		return int(x.max()+2)
	bin_width = 2*(q75 - q25)*len(x)**(-1/3)
	return int(round((x.max() - x.min())/bin_width))

fig, ax = plt.subplots(len(snakemake.input),1, sharex=True, figsize=(10,10))
#plt.subplots_adjust(hspace=0.5)
for i, filename in enumerate(snakemake.input):
	fn = filename.split("/")[-1].rsplit(".", 1)[0]
	print(fn)
	if fn == 'umi_per_gem.raw':
		df = pd.read_csv(snakemake.input[i], sep=r'\s+', header=None, names=['gem', 'freq'])
		title = 'Raw barcode reads'
		x = df.freq
		a,b = [],[]
		upper_limit = max(x)
	elif fn == 'umi':
		df = pd.read_csv(snakemake.input[i], sep='\t', usecols=['gem','A_N6','B_N6'])
		df['umi'] = df.A_N6 + df.B_N6

		x = df.groupby('gem').umi.size()
		title = 'Mapped barcodes'
		a,b = [],[]
	elif fn == 'mapping.clean.AKB.augmented':
		df = pd.read_csv(snakemake.input[i], usecols=['gem','template_id_mhc','umi_count_mhc','template_id_cd8','umi_count_cd8'])
		title = 'Credibly mapped barcodes'
		a = df.umi_count_mhc.dropna()
		b = df.umi_count_cd8.dropna()
		x = []
	elif fn == 'tcr_barcode.cleaned':
		df = pd.read_csv(snakemake.input[i], usecols=['gem','template_id_mhc','umi_count_mhc','template_id_cd8','umi_count_cd8'])
		title = 'Credibly mapped barcodes in GEMs also containing TCRs'
		a = df.umi_count_mhc.dropna()
		b = df.umi_count_cd8.dropna()
		x = []
	else:
		print('unkown filename')

	ax[i].hist(x, bins=get_bins(x), alpha=0.7, label='Any')
	ax[i].hist(b, bins=get_bins(b), alpha=0.7, label='TRB')
	ax[i].hist(a, bins=get_bins(a), alpha=0.7, label='TRA')
	ax[i].legend()
	ax[i].set_title(title, fontsize=10)
	ax[i].set_xscale('log')
	ax[i].set_xlim(0.91, upper_limit)

fig.text(0.5, 0.04, 'UMIs per GEM', ha='center')
fig.text(0.04, 0.5, 'Counts', va='center', rotation='vertical')
fig.suptitle('UMI distributions for GEMs containing barcodes')
#fig.set_xscale('log')
#fig.set_xlim(0.9, upper_limit)

plt.savefig(snakemake.output[0], bbox_inches='tight')
plt.show()
