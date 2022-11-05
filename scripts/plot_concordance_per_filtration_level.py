#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import re

from D_plot_specificity_matrix_utils import (filtering, calc_binding_concordance)

df = pd.read_csv(snakemake.input.table)
total_gems = df.shape[0]

i = 1
fig = plt.figure(figsize=(20,15))

for file_X, file_Y in zip(snakemake.input.file_X, snakemake.input.file_Y): #snakemake.input.conc_Z,  gems_Z , snakemake.input.gems_Z

    clonotype_fmt	 = file_X.split('/')[10]
    delta_umi		 = file_X.split('/')[11]
    filtration_level = file_X.split('/')[12]
    delta = int(re.search('(\d+)', delta_umi).group(1))
    
    print(i)
    print(clonotype_fmt, delta_umi, filtration_level)
    
    assert file_X.split('/')[11] == file_Y.split('/')[11]
    assert file_X.split('/')[12] == file_X.split('/')[12]

    X = np.loadtxt(file_X)
    Y = np.loadtxt(file_Y)
    #Zc = np.loadtxt(conc_Z)
    #Zg = np.loadtxt(gems_Z)
    #Zg = Zg / total_gems
    Zc = X * 0
    Zg = X * 0

    for j, (bumis, tumis) in enumerate(zip(X,Y)):
        for k, (bmin, tmin) in enumerate(zip(bumis, tumis)):

            tmp_df = filtering(df, bmin, tmin, clonotype_fmt=clonotype_fmt, filtration=filtration_level, umi_delta=delta)
            tmp_df = calc_binding_concordance(tmp_df, 'ct')
            
            #assert len(tmp_df) == Zg[j,k]

            Zc[j,k] = tmp_df.binding_concordance.mean()
            Zg[j,k] = len(tmp_df)
            
    Zg = Zg / total_gems

    ax = fig.add_subplot(3, 3, i, projection='3d')
    ax.plot_surface(X, Y, Zc, rstride=1, cstride=1, cmap=cm.coolwarm, label='Average concordance')
    ax.plot_surface(X, Y, Zg, rstride=1, cstride=1, cmap='winter', label='Pct. GEM')
    ax.set_xticks(X[0].astype(int))
    ax.set_xlabel('min. BC UMI count')
    ax.set_ylabel('min. TCR UMI count')
    ax.set_zlabel('Average concordance\nFraction retained GEMs')
    ind = np.unravel_index(np.argmax(Zc, axis=None), Zc.shape)
    ax.set_title('%s\n%s\n%d GEMs at max concordance of %.2f (%i,%i)' %(delta_umi, filtration_level, Zg[ind]*total_gems, Zc[ind], X[ind], Y[ind])) 



    i += 1
plt.savefig(snakemake.output[0], bbox_inches='tight')
