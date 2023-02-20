#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import os
import json
import argparse

def get_argparser():
    """Return an argparse argument parser."""
    parser = argparse.ArgumentParser(prog = 'Augment merged data',
                                     description = 'Augment merged data')
    add_arguments(parser)
    return parser

def add_arguments(parser):
    parser.add_argument('--hto', required=True, help='Filepath to output of Seurat HTO analysis')
    parser.add_argument('--gex', required=True, help='Filepath to output of gene expression analysis')
    parser.add_argument('--tcrdb', required=True, help='Filepath to download of IEDB, VDJdb and MIRA specificities')
    parser.add_argument('--clones', required=False, default=None, help='Filepath to data of TCR clones')
    parser.add_argument('--report', required=False, default=None, help='Filepath for report of summary statistics')
    parser.add_argument('--data', required=True, help='Filepath of merged TCR and barcode data')
    parser.add_argument('--output', required=True, help='Filepath of output data')
    
def calc_binding_concordance(df, clonotype_fmt):
    gems_per_specificity        = df.groupby([clonotype_fmt,'peptide_HLA']).gem.count().to_dict()
    df['gems_per_specificity']  = df.set_index([clonotype_fmt,'peptide_HLA']).index.map(gems_per_specificity)
    gems_per_spec_hla_match     = df[df.HLA_match == True].groupby([clonotype_fmt, 'peptide_HLA']).gem.count().to_dict()
    df['gems_per_spec_hla_match'] = df.set_index([clonotype_fmt,'peptide_HLA']).index.map(gems_per_spec_hla_match)
    gems_per_clonotype          = df.groupby([clonotype_fmt]).gem.count().to_dict()
    df['gems_per_clonotype']    = df[clonotype_fmt].map(gems_per_clonotype)
    df['binding_concordance']   = df.gems_per_specificity / df.gems_per_clonotype
    df['hla_concordance']       = df.gems_per_spec_hla_match / df.gems_per_specificity
    df['hla_concordance']       = df.hla_concordance.fillna(0)
    return df

def calc_relative_umi(grp):
    return grp.umi_count_mhc / grp.umi_count_mhc.quantile(0.9, interpolation='lower')

def check_tcrdb(credible_df):
    """
    Check if any of the clonotypes correspond to known clonotypes from database of TCR specificities.
    The CDR3s of IEDB may have had C removed from sequence, and thus may be a substring of the 10x CDR3 sequences.
    For alpha and beta chain separately, find the matching substring in TCRdb (if any).
    Map the TCRdb peptide annotation to our data.
    Check if the TCRdb mapped peptide corresponds to the 10x annotated peptide.
    FYI, even public databases have 'cross-reactive' TCR annotations...
    """

    # List CDR3s
    all_A3 = '|'.join(tcrdb.cdr3_a.dropna().unique())
    all_B3 = '|'.join(tcrdb.cdr3_b.dropna().unique())
    
    # Find matching substring in TCRdb (if any)
    credible_df['cdr3_TRA_substr'] = credible_df.cdr3_TRA.fillna('').str.findall(all_A3).apply(lambda x: x[0] if len(x)==1 else np.nan)
    credible_df['cdr3_TRB_substr'] = credible_df.cdr3_TRB.fillna('').str.findall(all_B3).apply(lambda x: x[0] if len(x)==1 else np.nan)
    
    # Map
    dct = tcrdb.groupby(['cdr3_a','cdr3_b']).peptide.apply(list)
    credible_df['tcrdb_pep'] = credible_df.set_index(['cdr3_TRA_substr','cdr3_TRB_substr']).index.map(dct)
    
    # Check
    credible_df['tcrdb_check'] = credible_df.apply(lambda row: any([row.peptide==pep for pep in row.tcrdb_pep]) if row.tcrdb_pep==row.tcrdb_pep else np.nan, axis=1)

    return credible_df


try:
    merged_annotations = snakemake.input.dat
    HTO = snakemake.input.hto
    GEX = snakemake.input.gex
    TCRdb = snakemake.input.vdj
    clone_sequencing = snakemake.params[0]
    output = snakemake.output[0]
    report = snakemake.output[1]
except:
    parser = get_argparser()
    args = parser.parse_args()
    
    merged_annotations = args.data
    HTO = args.hto
    GEX = args.gex
    TCRdb = args.tcrdb
    clone_sequencing = args.clones
    output = args.output
    report = args.report

# # Import input data
gex = pd.read_csv(GEX, header=None, names=['gem'])
hto = pd.read_csv(HTO, skiprows=1, header=None, names=['gem','seurat','umi_count_hto','feature_hto','hto_max_id','hto_sec_id',
                                                       'hto_margin','hto_classification','hto_global_class','hash_id'])
credible_df = pd.read_csv(merged_annotations, low_memory=False)
tcrdb = pd.read_csv(TCRdb)
if clone_sequencing:
    cs_df = pd.read_csv(clone_sequencing, sep='\t')
    cs_df.drop_duplicates(subset=['amino_acid'], inplace=True)

rprt = {'tcr_barcode_merged_gems': list(set(credible_df.gem))}

##############################################################################################################################
#                                                            Clean                                                           #
##############################################################################################################################
# Remove NaN, None.
assert credible_df.gem.isna().sum() == 0
credible_df.dropna(subset=['gem','ct','peptide_HLA'], inplace=True)
rprt['rm_nan_mhc_gems'] = list(set(credible_df.gem))

# Remove GEMs annotated with epitope 0
credible_df.drop(credible_df[credible_df.epitope == '0'].index, inplace=True)
rprt['no-0-epitope_gems'] = list(set(credible_df.gem))

###############################################################################################################################
#                                        Augment (Relative UMI count within clonotype)                                        #
###############################################################################################################################
for ct, grp in credible_df.groupby('ct'):
    umi_rel = grp.umi_count_mhc / grp.umi_count_mhc.quantile(0.9, interpolation='lower') #.max()
    credible_df.loc[umi_rel.index, 'umi_count_mhc_rel'] = umi_rel

###############################################################################################################################
#                                                     Check clones & Merge                                                    #
###############################################################################################################################
if clone_sequencing:
    credible_df = pd.merge(credible_df, cs_df, how='left', left_on='cdr3_TRB', right_on='amino_acid')

credible_df = calc_binding_concordance(credible_df, 'ct')

###############################################################################################################################
#                                                         Sanity check                                                        #
###############################################################################################################################
# Only include GEMs where there is a single clonotype 
if not credible_df.groupby(['gem']).ct.nunique().eq(1).all():
    print(credible_df.groupby(['gem']).ct.nunique().eq(1))


###############################################################################################################################
#                                                         Check TCRdb                                                         #
###############################################################################################################################
credible_df = check_tcrdb(credible_df)

###############################################################################################################################
#                                                       ADD HTO analysis                                                      #
###############################################################################################################################
df = pd.merge(credible_df, hto, how='left', on='gem')

###############################################################################################################################
#                                                       ADD GEX analysis                                                      #
###############################################################################################################################
df['gex'] = df.gem.isin(gex.gem)

###############################################################################################################################
#                                                            Write                                                            #
###############################################################################################################################
df.to_csv(output, index=False)

if report:
    with open(report, 'w') as fh:
        json.dump(rprt, fh)





