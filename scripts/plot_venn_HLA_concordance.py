#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn2

fig, axes = plt.subplots(1,len(snakemake.input))
for filename, ax in zip(snakemake.input, axes):
	# # Load
	df = pd.read_csv(filename)

	# # Split data
	gem_mhc = df[~df.template_id_mhc.isna()].gem
	gem_cd8_t = df[df.HLA_match == True].gem
	gem_cd8_f = df[df.HLA_match == False].gem

	# # Plot
	a, b, c = set(gem_mhc), set(gem_cd8_t), set(gem_cd8_f)

	ab = a.intersection(b)
	ac = a.intersection(c)
	bc = b.intersection(c)

	abc = len(ab.intersection(c))
	ab_c = len(ab) - abc
	ac_b = len(ac) - abc
	bc_a = len(bc) - abc
	a_bc = len(a) - ab_c - ac_b - abc
	b_ac = len(b) - ab_c - bc_a - abc
	c_ab = len(c) - ac_b - bc_a - abc

	venn3(subsets = (a_bc, b_ac, ab_c, c_ab, ac_b, bc_a, abc),
		  set_labels = ('MHC', 'Ab (HLA match)', 'Ab (HLA mismatch)'), ax=ax)
	if 'imput' in filename:
		ax.set_title('After imputations')
	elif 'clean' in filename:
		ax.set_title('After cleaning\nBefore imputations')
	else:
		ax.set_title('Before cleaning')
title = "GEMs containing MHC and Ab BCs in HLA concordance"
fig.suptitle(title)
fig.savefig(snakemake.output[0], bbox_inches='tight')

