#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
from ast import literal_eval
import re

def HLA_cd8_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace(",", "").replace("'","").split(" ")

def cdr3_lst_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace("'","").split(" ")

def epitope_converter(x):
    #define format of datetime
    return [y for y in x.replace("[","").replace("]","").replace("\n","").split("'") if (y != '') & (y != ' ')]

def peptide_hla_converter(x):
    return re.findall("\w+\s{1}\w{1}\d+", x.replace("[","").replace("]","").replace("\n","").replace("'",""))

converters = {'peptide_HLA_lst': peptide_hla_converter,'epitope_lst': epitope_converter,
              'umi_count_lst_mhc': literal_eval,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'HLA_lst_mhc': cdr3_lst_converter,'HLA_cd8': HLA_cd8_converter} #

# # Load
df = pd.read_csv(snakemake.input[0], converters=converters)

selected_clonotypes = df.groupby('ct').size().sort_values(ascending=False).head(10)
all_peptide_HLA = df.peptide_HLA_lst.explode()

for ct in selected_clonotypes.index:
    mat = pd.DataFrame(index=all_peptide_HLA.unique(), columns=df.gem.unique())
    df[df.ct == ct].apply(lambda row: sum_umi(row), axis=1)
    summed_umis = mat.sum(axis=1).sort_values().sort_values(ascending=False).head(10)
    summed_gems = (mat > 0).sum(axis=1)
    annotated_gems = df[df.ct == ct].groupby('peptide_HLA').size()
    all_gems_hla = pd.concat([annotated_gems,
                              pd.Series(0, index=all_peptide_HLA[~all_peptide_HLA.isin(annotated_gems.index)].unique())])
    
    plt.figure(figsize=(4,3)) #, bins=mat.T[peptide_HLA].dropna().max()
    for peptide_HLA in summed_umis.index:
        plt.hist(mat.T[peptide_HLA].dropna(), alpha=0.2, label='%s %s %d (%d) (%d)' %(peptide_HLA.split(' ')[1], peptide_HLA.split(' ')[0], summed_umis[peptide_HLA], summed_gems[peptide_HLA], all_gems_hla[peptide_HLA]))
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', title='HLA peptide UMI-sum (found in #GEMs) (ann. in #GEMs)')
    plt.xlabel('UMI count')
    plt.ylabel('GEM count')
    plt.title('Clonotype %d (%d GEMs)\n%s - %s\n%s - %s' %(ct, selected_clonotypes[ct],
                                               df[df.ct == ct].genes_TRA.unique()[0], df[df.ct == ct].cdr3_TRA.unique()[0],
                                               df[df.ct == ct].genes_TRB.unique()[0], df[df.ct == ct].cdr3_TRB.unique()[0]))

    r = re.compile('.*ct_%.pdf'%ct)
    output = list(filter(r.match, snakemake.output))
    assert len(output) == 1
    plt.savefig(output[0], bbox_to_anchor='tight')
    
    print(df[df.ct == ct].cdr3_TRA.unique())
    print(df[df.ct == ct].genes_TRA.unique())
    print(df[df.ct == ct].cdr3_TRB.unique())
    print(df[df.ct == ct].genes_TRB.unique())
