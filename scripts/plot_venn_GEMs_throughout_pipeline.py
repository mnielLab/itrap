#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import json
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn2

# OBS OBS OBS!
# Scale diagrams to size
# https://stackoverflow.com/questions/34183870/matplotlib-venn-multiple-venn-diagrams-same-scale

def add_venn(index,a,b, a_key, b_key, rows, cols):
	ab = len(a.intersection(b))
	a_b = len(a) - ab
	b_a = len(b) - ab

	#print('figure number: %d, %d, %d' %(i,j, i+j+1+(i*2)))
	fig.add_subplot(rows,cols,index)
	venn2(subsets = (a_b, b_a, ab), set_labels = ('TCR', 'BC'))
	title = '%d GEMs, %d %s (TCR), %d %s (BC)' %(len(a.union(b)), len(a), a_key, len(b), b_key)
	plt.title(title)

def split_label(label):
	label = label.replace('total_', '')
	la = label.split('_')
	if len(la) < 2:
		return la[0]
	if len(la) > 4:
		l = ' '.join(la).split(' ', 2)
		print(l)
		return '%s\n%s' %(' '.join(l[:2]), l[-1])
	l = ' '.join(la).split(' ', 1)
	return '%s\n%s' %(l[0], l[-1])

TCR_GEMS = snakemake.input.TCR_GEMS
BAR_GEMS = snakemake.input.BAR_GEMS

tcr_dct = dict()
bar_dct = dict()

# Load TCR GEMS
for tcr_file in TCR_GEMS:
	tcr_f = tcr_file.split('.')
	if tcr_f[-1] == 'json':
		with open(tcr_file) as json_file:
			tcr_dct.update(json.load(json_file))
	else:
		tcr_v = np.loadtxt(tcr_file, dtype='U32')
		tcr_k = tcr_f[-2]
		tcr_dct[tcr_k] = tcr_v

# Load Barcode GEMs
for bar_file in BAR_GEMS:
	bar_f = bar_file.split('.')
	if bar_f[-1] == 'json':
		with open(bar_file) as json_file:
			bar_dct.update(json.load(json_file))
		del bar_dct['total_umi_mapped_gems']
		del bar_dct['full_length_umi_gems']
	else:
		bar_v = np.loadtxt(bar_file, dtype='U32')
		bar_k = bar_f[-2]
		bar_dct[bar_k] = bar_v

rows = len(tcr_dct)
cols = len(bar_dct)
#fig, ax = plt.subplots(rows, cols, figsize=(15,10))
fig = plt.figure(figsize=(15,15))
index = 0
for i, (tcr_k, tcr_v) in enumerate(tcr_dct.items()):
	for j, (bar_k, bar_v) in enumerate(bar_dct.items()):
		index += 1

		#print(bar_k)

		# Get Venn
		a, b = set(tcr_v), set(bar_v)
		ab = len(a.intersection(b))
		a_b = len(a) - ab
		b_a = len(b) - ab
		#add_venn(index, set(tcr_v), set(bar_v), tcr_k, bar_k, rows, cols)

		# Plot
		fig.add_subplot(rows,cols,index)
		if i == 0:
			plt.text(0, 1.2, '%d %s' %(len(b), split_label(bar_k)), ha='center', va='bottom', rotation='horizontal') #, transform=ax.transAxes
		if j == 0:
			plt.text(-1, 0, '%d %s' %(len(a), ' '.join(tcr_k.split('_'))), ha='right', va='center', rotation='vertical') #, transform=ax.transAxes
		v = venn2(subsets = (a_b, b_a, ab), set_labels = ('%d TCR' %a_b, '%d BC' %b_a))
		# setting the font size
		for t in v.set_labels:
			t.set_fontsize(6)
		# Remove numbers except the interception
		v.get_label_by_id('100').set_text('')
		v.get_label_by_id('010').set_text('')
		for t in v.subset_labels: # number font size
			if t is not None:
				t.set_fontsize(10)
		title = '%d GEMs' %len(a.union(b))
		plt.title(title)
		
plt.suptitle("GEMs containing TCR and/or BCs of MHC + Ab")
plt.savefig(snakemake.output[0], bbox_inches='tight')