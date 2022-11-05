#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import re
from ast import literal_eval
import itertools

from D_plot_specificity_matrix_utils import SpecificityMatrix, calc_binding_concordance

# # Args
def cdr3_lst_converter(x):
    #define format of datetime
    return x.replace("[","").replace("]","").replace("'","").split(" ")

def peptide_hla_converter(x):
    return re.findall("\w+\s{1}\w{1}\d+", x.replace("[","").replace("]","").replace("'",""))

def literal_converter(val):
    try:
        return literal_eval(val)
    except SyntaxError:
        return np.nan
    except ValueError:
        return np.nan

converters = {'peptide_HLA_lst': peptide_hla_converter,
              'cdr3_lst_TRA': cdr3_lst_converter,
              'cdr3_lst_TRB': cdr3_lst_converter,
              'umi_count_lst_mhc': literal_eval,
              'umi_count_lst_TRA': literal_converter,
              'umi_count_lst_TRB': literal_converter}

class Impute():
    concordance_threshold = 0.5
    threshold_delta_umi = 4 #100 with new delta UMI computation the minimum delta value for a singlet is 4=1/0.25
    colR = ['gem', 'ct', 'cdr3_TRA','umi_count_TRA', 'cdr3_TRB','umi_count_TRB', 'peptide_HLA', 'umi_count_mhc', 'peptide', 'HLA_cd8', 'umi_count_cd8']
    cols = ['gem', 'cdr3_TRA', 'umi_count_TRA', 'cdr3_TRB', 'umi_count_TRB', 'peptide', 'umi_count_mhc', 'peptide_HLA', 'ct', 'epitope']

    def __init__(self, df):
        self.df = df
        self.stats = pd.DataFrame(index=['potential', 'imputations', 'percent', 'corrected'])

    @staticmethod
    def prepare_dataframe(df):
        df = calc_binding_concordance(df, 'ct')
        df['chain_count_TRA'] = df.apply(lambda x: len(x.cdr3_lst_TRA) if x.cdr3_lst_TRA[0] != '' else 0, axis=1)
        df['chain_count_TRB'] = df.apply(lambda x: len(x.cdr3_lst_TRB) if x.cdr3_lst_TRB[0] != '' else 0, axis=1)
        return df

    def compute_stats(self, query, query_filter, hit):
        if query == 'TRA_NA':
            variable = 'cdr3_TRA'
        elif query == 'TRB_NA':
            variable = 'cdr3_TRB'
        elif query == 'TRA':
            variable = 'cdr3_TRA'
        elif query == 'TRB':
            variable = 'cdr3_TRB'
        elif query == 'PEP':
            variable = 'peptide_HLA'
        elif query == 'ref':
            variable = 'peptide_HLA'

        self.stats[query] = [sum(query_filter),
                             len(hit),
                             "%.2f" %round(len(hit)/sum(query_filter)*100, 2),
                             sum(self.df.loc[hit.index, variable] != hit.apply(pd.Series)[0])]

    def get_reference(self):
        """
        Query filter: we still dont have a great overlap between pMHC and Ab GEMs,
        so for the reference deselect GEMs where we know there is a mismatch.
        """
        query_filter = ((self.df.HLA_match != False) &
                        (self.df.single_barcode_cd8 != False) &
                        (self.df.multiplets_mhc == 1) &
                        (self.df.chain_count_TRA == 1) &
                        (self.df.chain_count_TRB == 1) &
                        (self.df.binding_concordance >= self.concordance_threshold))
        self.df.loc[query_filter, 'imputed_TRA'] = self.df[query_filter].cdr3_TRA
        self.df.loc[query_filter, 'imputed_TRB'] = self.df[query_filter].cdr3_TRB
        self.df.loc[query_filter, 'imputed_PEP'] = self.df[query_filter].peptide_HLA
        self.df.loc[query_filter, 'imputed_CDX'] = self.df[query_filter].HLA_cd8
        self.df.loc[query_filter, 'imputed_umi_TRA'] = self.df[query_filter].umi_count_TRA
        self.df.loc[query_filter, 'imputed_umi_TRB'] = self.df[query_filter].umi_count_TRB
        self.df.loc[query_filter, 'imputed_umi_PEP'] = self.df[query_filter].umi_count_mhc
        self.df.loc[query_filter, 'imputation_category'] = 'unique'
        self.reference = self.df[query_filter].drop_duplicates(subset=['cdr3_TRA', 'cdr3_TRB', 'peptide_HLA']).loc[:, self.colR]

        self.compute_stats('ref', query_filter, hit=self.df[query_filter].peptide_HLA)

    def get_hit(self, query):
        annotate = self.get_annotation_fct(query)
        query_filter = self.get_filter(query)

        hit = self.df[query_filter].apply(lambda row: annotate(row), axis=1).dropna() #self.threshold_delta_umi
        print(query)
        print(hit.empty)
        print(query_filter)

        if not hit.empty:
            q = query[:3]
            self.df[['imputed_{}'.format(q), 'imputed_umi_{}'.format(q)]] = hit.apply(pd.Series)
            self.df.loc[hit.index, 'imputation_category'] = query

            self.compute_stats(query, query_filter, hit)

    def get_filter(self, query):
        if query == 'TRA_NA':
            return ((self.df.multiplets_mhc == 1) & (self.df.chain_count_TRA == 0) & (self.df.chain_count_TRB == 1))
        elif query == 'TRB_NA':
            return ((self.df.multiplets_mhc == 1) & (self.df.chain_count_TRA == 1) & (self.df.chain_count_TRB == 0))
        elif query == 'TRA':
            return ((self.df.multiplets_mhc == 1) & (self.df.chain_count_TRA > 1) & (self.df.chain_count_TRB == 1))
        elif query == 'TRB':
            return ((self.df.multiplets_mhc == 1) & (self.df.chain_count_TRA == 1) & (self.df.chain_count_TRB > 1))
        elif query == 'PEP':
            return ((self.df.multiplets_mhc > 1) & (self.df.chain_count_TRA == 1) & (self.df.chain_count_TRB == 1))

    def get_annotation_fct(self, query):
        if query == 'TRA_NA':
            return self.annotate_TRA_NA
        elif query == 'TRB_NA':
            return self.annotate_TRB_NA
        elif query == 'TRA':
            return self.annotate_TRA
        elif query == 'TRB':
            return self.annotate_TRB
        elif query == 'PEP':
            return self.annotate_PEP
        elif query == 'CDX':
            return self.annotate_CDX

    def annotate_TRA_NA(self, row):
        assert row.cdr3_TRA != np.nan
        hit = self.reference[(self.reference.peptide_HLA == row.peptide_HLA) & (self.reference['cdr3_TRB'] == row['cdr3_TRB'])].copy()
        if len(hit) == 1:
            return hit[['cdr3_TRA', 'umi_count_TRA']].values[0]
        return None

    def annotate_TRB_NA(self, row):
        assert row.cdr3_TRB != np.nan
        hit = self.reference[(self.reference.peptide_HLA == row.peptide_HLA) & (self.reference['cdr3_TRA'] == row['cdr3_TRA'])].copy()
        if len(hit) == 1:
            return hit[['cdr3_TRB', 'umi_count_TRB']].values[0]
        return None

    def annotate_TRA(self, row):
        query = self.reference[(self.reference.peptide_HLA == row.peptide_HLA) & (self.reference['cdr3_TRB'] == row['cdr3_TRB'])].copy()
        
        if row.delta_umi_TRA < self.threshold_delta_umi:
            query['hit'] = query['cdr3_TRA'].isin(row['cdr3_lst_TRA'])
        else:
            return None
        
        hit = query[query.hit].copy()

        if len(hit) == 1:
            element = hit['cdr3_TRA'].values[0]
            index = row['cdr3_lst_TRA'].index(element)
            return element, row['umi_count_lst_TRA'][index]
        elif len(hit) > 1:
            return row['cdr3_TRA'], row['umi_count_TRA']       
        return None

    def annotate_TRB(self, row):
        query = self.reference[(self.reference.peptide_HLA == row.peptide_HLA) & (self.reference['cdr3_TRA'] == row['cdr3_TRA'])].copy()
        
        if row.delta_umi_TRB < self.threshold_delta_umi:
            query['hit'] = query['cdr3_TRB'].isin(row['cdr3_lst_TRB'])
        else:
            return None
        
        hit = query[query.hit].copy()

        if len(hit) == 1:
            element = hit['cdr3_TRB'].values[0]
            index = row['cdr3_lst_TRB'].index(element)
            return element, row['umi_count_lst_TRB'][index]
        elif len(hit) > 1:
            return row['cdr3_TRB'], row['umi_count_TRB']
        return None

    def annotate_PEP(self, row):
        query = self.reference[(self.reference.cdr3_TRA == row.cdr3_TRA) & (self.reference.cdr3_TRB == row.cdr3_TRB)].copy()
        
        if row.delta_umi_mhc < self.threshold_delta_umi:
            query['hit'] = query['peptide_HLA'].isin(row['peptide_HLA_lst'])
        else:
            return None
        
        hit = query[query.hit].copy()

        if len(hit) == 1:
            element = hit['peptide_HLA'].values[0]
            index = row['peptide_HLA_lst'].index(element)
            return element, row['umi_count_lst_mhc'][index]
        elif len(hit) > 1:
            print(hit)
            return row['peptide_HLA'], row['umi_count_mhc']
        return None

    def annotate_CDX(self):
        """ 
        This is not done by row, but for entire df at once (unlike previous functions).
        When imputing CDX (sample HLA) it doesn't matter which CDR3 we have. We just need a reliable relation between pMHC and sample HLAs.
        Function handles two cases:
        1) when HLA match is True: check if any other sample ID fits in?
        2) when HLA match is False: check if any oterh sample ID fits in!
        """
        sample = self.df[self.df.single_barcode_cd8 == False]
        query = self.reference[(self.reference.peptide_HLA == row.peptide_HLA)].copy()
        
        if row.delta_umi_mhc < self.threshold_delta_umi:
            query['hit'] = query['peptide_HLA'].isin(row['peptide_HLA_lst'])
        else:
            return None
        
        hit = query[query.hit].copy()

        if len(hit) == 1:
            element = hit['peptide_HLA'].values[0]
            index = row['peptide_HLA_lst'].index(element)
            return element, row['umi_count_lst_mhc'][index]
        elif len(hit) > 1:
            print(hit)
            return row['peptide_HLA'], row['umi_count_mhc']
        return None

    def modify(self):
        df = self.df.copy()

        self.df['cdr3_TRA'] = np.where(df.imputed_TRA.isna(), df.cdr3_TRA, df.imputed_TRA)
        self.df['cdr3_TRB'] = np.where(df.imputed_TRB.isna(), df.cdr3_TRB, df.imputed_TRB)
        self.df['peptide_HLA'] = np.where(df.imputed_PEP.isna(), df.peptide_HLA, df.imputed_PEP)
        self.df['umi_count_TRA'] = np.where(df.imputed_umi_TRA.isna(), df.umi_count_TRA, df.imputed_umi_TRA)
        self.df['umi_count_TRB'] = np.where(df.imputed_umi_TRB.isna(), df.umi_count_TRB, df.imputed_umi_TRB)
        self.df['umi_count_mhc'] = np.where(df.imputed_umi_PEP.isna(), df.umi_count_mhc, df.imputed_umi_PEP)

        self.df['cdr3_comb'] = self.df.cdr3_TRA.fillna('') + self.df.cdr3_TRB.fillna('')

        self.df['ct'] = self.assign_clonotype(df)

    def assign_clonotype(self, df):
        clonotype_variables = ['cdr3_TRA','cdr3_TRB']
        df.loc[:, clonotype_variables] = df.loc[:, clonotype_variables].fillna('unknown')
        new_clonotype = df.groupby(clonotype_variables).gem.unique().to_frame()
        new_clonotype['n_gems'] = new_clonotype.gem.apply(len)
        new_clonotype.sort_values(by='n_gems', ascending=False, inplace=True)
        dct = new_clonotype.to_dict()['gem']
        for i, k in enumerate(dct.keys(), start=1): 
            dct[k] = i
        return df.set_index(clonotype_variables).index.map(dct)
    
# ## Data
if __name__ == '__main__':
    INPUT = snakemake.input.df
    OUTPUT = snakemake.output.output
    REPORT = snakemake.output.report

    # # Import
    df = pd.read_csv(INPUT, converters=converters, low_memory=False)
    df = Impute.prepare_dataframe(df)

    inst = Impute(df)
    inst.get_reference() # reference should also include Ab HLAs
    # impute Ab HLA for GEMs with single MHC but multiple sample IDs (GEMs that would have been part of reference if I didn't filter on HLA-match)
    # update reference 
    inst.get_hit('TRA_NA')
    inst.get_hit('TRB_NA')
    inst.get_hit('TRA')
    inst.get_hit('TRB')
    inst.get_hit('PEP')
    inst.modify()

    # ## Write to excel
    inst.df.to_csv(OUTPUT, index=False)
    inst.stats.to_csv(REPORT, index=False)


