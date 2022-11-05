#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from itertools import groupby

#def group_files(snakemake_input):
#	for key, group in groupby(snakemake_input, lambda x: x.split('/')[10]):
#		#files = list()
#		files = dict()
#		for k, g in groupby(list(group), lambda x: x.split('/')[-1]):
#			#files.append(list(g))
#			files[k] = list(g)
#		yield zip(files['X'], files['Y'], files['Z_siml'], files['Z_gems']) # without list?

df = pd.read_csv(snakemake.input.table)
total_gems = df.shape[0]

#for fmt_files in group_files(snakemake.input):
i = 1
fig = plt.figure(figsize=(20,15))
for file_X, file_Y, siml_Z, gems_Z in zip(snakemake.input.file_X, snakemake.input.file_Y, snakemake.input.siml_Z, snakemake.input.gems_Z): #fmt_files:

	clonotype_fmt	 = file_X.split('/')[10]
	delta_umi		 = file_X.split('/')[11]
	filtration_level = file_X.split('/')[12]
	print(i)
	print(clonotype_fmt, delta_umi, filtration_level)
	
	assert file_X.split('/')[11] == file_Y.split('/')[11]
	assert file_X.split('/')[12] == file_X.split('/')[12]

	X = np.loadtxt(file_X)
	Y = np.loadtxt(file_Y)
	Zc = np.loadtxt(siml_Z)
	Zg = np.loadtxt(gems_Z)
	Zg = Zg / total_gems

	ax = fig.add_subplot(3, 3, i, projection='3d')
	ax.plot_surface(X, Y, Zg, rstride=1, cstride=1, cmap='winter', alpha=0.5, label='Pct. GEM')
	ax.plot_surface(X, Y, Zc, rstride=1, cstride=1, cmap=cm.coolwarm, edgecolor='black', label='Average concordance')
	ax.set_xticks(X[0].astype(int))
	ax.set_xlabel('min. BC UMI count')
	ax.set_ylabel('min. TCR UMI count')
	ax.set_zlabel('Fraction of significant groups\nFraction retained GEMs', labelpad=10)
	ind = list(zip(*np.where(Zc == Zc.max())))[-1] #np.unravel_index(np.argmax(Zc, axis=None), Zc.shape)
	ax.set_title('%s\n%s\n%d GEMs at sign. frac. %.2f (%i,%i)' %(delta_umi, filtration_level, Zg[ind]*total_gems, Zc[ind], X[ind], Y[ind])) 

	i += 1

#plt.subplots_adjust(left= 0.001, wspace=0.003) #, right=0.012
plt.tight_layout(w_pad=0.5) #pad=0.4, 
plt.savefig(snakemake.output[0], bbox_inches='tight')
