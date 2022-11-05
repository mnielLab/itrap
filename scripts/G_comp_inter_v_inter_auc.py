#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from ast import literal_eval
import re
from scipy import stats
from random import sample
#import scipy.stats as stats
import random
import re
import statistics
from scipy.stats import wilcoxon
from sklearn.metrics import roc_curve, auc


###########################################################
#                        Functions                        #
###########################################################

def notnan(x):
    return x == x

##########################################################
#                        IN/Output                       #
##########################################################
SIM_INPUT = snakemake.input.sim
AUC_OUTPUT = snakemake.output.auc
#filter_set = snakemake.params.flt

##########################################################
#                          Load                          #
##########################################################
test_df = pd.read_csv(SIM_INPUT, header=None, names=['plateau', 'rnd_sample','filtering','peptide','ct','score','cdr3s'])

##########################################################
#               Compute AUC for each sample              #
##########################################################
test_df['label'] = np.where(test_df.plateau == 'intra', 2, 0)

for s, s_grp in test_df.dropna().groupby('rnd_sample'): # dropna should be superflouos 
    for f, f_grp in s_grp.groupby('filtering'):
        fpr, tpr, _ = roc_curve(f_grp.label, f_grp.score, pos_label=2)
        test_df.loc[f_grp.index, 'AUC'] = auc(fpr, tpr)
        test_df.loc[f_grp.index, 'AUC 0.1'] = auc(fpr[fpr < 0.1], tpr[fpr < 0.1])/ 0.1

auc_df = test_df.dropna().melt(id_vars=['rnd_sample', 'filtering'], value_vars=['AUC','AUC 0.1'])
auc_df['palette'] = auc_df.filtering.map(dict(zip(labels,palette)))

##########################################################
#                  Write output to file                  #
##########################################################
#test_df.to_csv(SIM_OUTPUT, index=False, header=False)
auc_df.to_csv(AUC_OUTPUT, index=False)



