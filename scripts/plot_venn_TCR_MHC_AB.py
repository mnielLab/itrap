#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_unweighted
from matplotlib_venn import venn2

def make_plain(v, weighted=True):
    """
    Convert venn to black and white
    """
    idxs = ['100', '010', '110', '001', '101', '011', '111']
    
    for i, s in zip(idxs, sets):
        try:
            v.get_patch_by_id(i).set_color('none')
            v.get_patch_by_id(i).set_edgecolor('black')
            v.get_patch_by_id(i).set_alpha(1)
            v.get_patch_by_id(i).set_lw(1.5)
        except:
            print(i, s)
            
    #v.get_patch_by_id('100').set_color('none')
    #v.get_patch_by_id('010').set_color('none')
    #v.get_patch_by_id('110').set_color('none')
    #v.get_patch_by_id('101').set_color('none')
    #v.get_patch_by_id('011').set_color('none')
    #v.get_patch_by_id('001').set_color('none')
    #v.get_patch_by_id('111').set_color('none')
    #
    #v.get_patch_by_id('100').set_edgecolor('black')
    #v.get_patch_by_id('010').set_edgecolor('black')
    #v.get_patch_by_id('110').set_edgecolor('black')
    #v.get_patch_by_id('101').set_edgecolor('black')
    #v.get_patch_by_id('011').set_edgecolor('black')
    #v.get_patch_by_id('001').set_edgecolor('black')
    #v.get_patch_by_id('111').set_edgecolor('black')
    #v.get_patch_by_id('111').set_edgecolor('black')
    #
    #v.get_patch_by_id('100').set_alpha(1)
    #v.get_patch_by_id('010').set_alpha(1)
    #v.get_patch_by_id('110').set_alpha(1)
    #v.get_patch_by_id('101').set_alpha(1)
    #v.get_patch_by_id('011').set_alpha(1)
    #v.get_patch_by_id('001').set_alpha(1)
    #v.get_patch_by_id('111').set_alpha(1)
    #
    #v.get_patch_by_id('100').set_lw(1.5)
    #v.get_patch_by_id('010').set_lw(1.5)
    #v.get_patch_by_id('110').set_lw(1.5)
    #v.get_patch_by_id('101').set_lw(1.5)
    #v.get_patch_by_id('011').set_lw(1.5)
    #v.get_patch_by_id('001').set_lw(1.5)
    #v.get_patch_by_id('111').set_lw(1.5)
    
    for text in v.set_labels:
        text.set_fontsize(16)
        
    for text in v.subset_labels:
        if text is not None:
            text.set_fontsize(14)

    if weighted == False:
        #Move the numbers in the circles  
        pos = v.get_label_by_id("101").get_position()
        v.get_label_by_id("101").set_position((pos[0], -0.23))
        pos = v.get_label_by_id("011").get_position()
        v.get_label_by_id("011").set_position((pos[0], -0.23))

# # Load
df = pd.read_csv(snakemake.input[0])
tcr_df = pd.read_csv(snakemake.input[1])

label = snakemake.params.label

if 'tot' not in label:
    l = 1 if 'pos' in label else 0
    tcr_df = tcr_df[tcr_df.label == l].copy()

# # Split data
gem_tcr = tcr_df.gem
gem_mhc = df[df.umi_count_mhc.notnull()].gem #df.template_id_mhc.notnull()
gem_cd8 = df[df.umi_count_cd8.notnull()].gem #df.template_id_cd8.notnull()

# # Plot
a, b, c = set(gem_tcr), set(gem_mhc), set(gem_cd8)

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
sets = (a_bc, b_ac, ab_c, c_ab, ac_b, bc_a, abc) # 100, 010, 110, 001, 101, 011, 111

v = venn3(subsets = sets, set_labels = ('TCR', 'pMHC', 'Hashing Ab'))
make_plain(v, weighted=True)
#title = "GEMs containing TCR, pMHC, and/or Ab BCs"
#plt.title(title)
plt.savefig(snakemake.output.w, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

v = venn3_unweighted(subsets = sets, set_labels = ('TCR', 'pMHC', 'Hashing Ab'))
make_plain(v, weighted=False)
#title = "GEMs containing TCR, pMHC, and/or Ab BCs"
#plt.title(title)
plt.savefig(snakemake.output.u, bbox_inches='tight', dpi=300)
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window