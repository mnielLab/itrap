#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import json
import itertools
import argparse

def get_argparser():
    """Return an argparse argument parser."""
    parser = argparse.ArgumentParser(prog = 'Clean & augment TCR data',
                                     description = 'Reformats Cellranger TCR output and re-annotates clonotypes')
    add_arguments(parser)
    return parser

def add_arguments(parser):
    parser.add_argument('--contig', required=True,
                        help='Cellranger output typically found in cellranger/outs/multi/vdj_t/all_contig_annotations.csv')
    parser.add_argument('--consensus', required=True,
                        help=('Cellranger output typically found in '
                              'cellranger/outs/per_sample_outs/cellranger_{total}/vdj_t/consensus_annotations.csv'))
    parser.add_argument('--output', required=True, help='output filename')
    parser.add_argument('--report', required=False, default=None, help='filename for reporting script statistics')


# # Args
keep_only_is_cells = False
keep_only_high_confidence = False
keep_only_full_length = True
keep_only_productive = True
keep_only_unamibiguous_gems = False

# # Load arguments
try:
    CONTIG = snakemake.input.contig
    CLONOTYPES = snakemake.params.clonot
    OUTPUT = snakemake.output.output
    report = snakemake.output.report
except: #NameError
    parser = get_argparser()
    args = parser.parse_args()
    CONTIG = args.contig
    CLONOTYPES = args.consensus
    OUTPUT = args.output
    report = args.report


###########################################
#                  Load                   #
###########################################
"""
Load data from positive an negative stainings.
If only positive stainings
"""
df = pd.read_csv(CONTIG)

# clone_df contains the consensus annotations for genes, sequence, cdrs, positions... 
# each clonotype has 2 rows in the df: TRA & TRB
clone_df = pd.read_csv(CLONOTYPES, header=0)

###########################################
#               Preprocess                #
###########################################
# ## Rename
df.rename(columns={"barcode" : "gem"}, inplace=True)
df.rename(columns={"raw_clonotype_id" : "clonotype"}, inplace=True)
df.rename(columns={"umis" : "umi_count"}, inplace=True)
df.rename(columns={"reads" : "read_count"}, inplace=True)

df.clonotype.replace('None',np.nan)
df.clonotype = np.where(df.clonotype == None, np.nan, df.clonotype)

rprt = {'total_tcr_gems': list(set(df.gem))}

# ## Filter data
flt_report = 'filter'
if keep_only_is_cells:
    print('keep_only_is_cells')
    df = df[df.is_cell == True]
    flt_report += '_is_cell'
if keep_only_high_confidence:
    print('keep_only_high_confidence')
    df = df[df.high_confidence == True]
    flt_report += '_high_confidence'
if keep_only_full_length:
    print('keep_only_full_length')
    df = df[df.full_length == True]
    flt_report += '_full_length'
if keep_only_productive:
    print('keep only productive')
    df = df[df.productive == True]
    #df = df[df.productive == 'True']
    flt_report += '_productive'
if keep_only_unamibiguous_gems:
    print('Keep only unambiguous gems')
    df = df.groupby(['gem', 'chain']).filter(lambda x: len(x) == 1)
    flt_report += '_unambiguous'

rprt[flt_report] = list(set(df.gem))

###########################################
#                  Clean                  #
###########################################
# Remove GEMs containing unnatural CDR3 sequences
expected_cdr3 = df.cdr3.str.contains('^[GPAVLIMCFYWHKRQNEDST]+$', regex=True, na=False)
df = df[expected_cdr3]

###########################################
#                  Genes                  #
###########################################
# Define a clone by: CDR3, V-gene & J-gene, sep=;
df['genes'] = df.replace([None], ['']).fillna('').apply(lambda x: ';'.join(x[['cdr3','v_gene','j_gene']]), axis=1)
clone_df['genes'] = clone_df.replace([None], ['']).fillna('').apply(lambda x: ';'.join(x[['cdr3','v_gene','j_gene']]), axis=1)

###########################################
#             Clone Libraries             #
###########################################
# 1. Dict of all duplicated clones --> replace duplicates
# 2. Dict of unique clonotype labels --> impute unannotated GEMs

## Find clonotype duplicates and unique dict of genes per clonotype label
clone_a = clone_df[clone_df.chain == 'TRA'].copy()
clone_b = clone_df[clone_df.chain == 'TRB'].copy()

# df index: clonotypes, 2 columns: genes_TRA & genes_TRB
clone1 = pd.merge(clone_a.groupby('clonotype_id').genes.unique().to_frame(),
                   clone_b.groupby('clonotype_id').genes.unique().to_frame(),
                   how='outer', left_index=True, right_index=True, suffixes=['_TRA','_TRB'])
# split clonotypes with multiple genes within one annotation
clone1 = clone1.explode('genes_TRA').explode('genes_TRB')
"""
                genes_TRA                               genes_TRB
clonotype_id
clonotype1      CAASQNEKLTF;TRAV21;TRAJ48               CASRIGAAGNSPLHF;TRBV27;TRBJ1-6
clonotype10     CAASEPAKIPKAAGNKLTF;TRAV29/DV5;TRAJ17   CSVEVPGKVFRTEAFF;TRBV29-1;TRBJ1-1
clonotype10     CAVETSGSRLTF;TRAV39;TRAJ58              CSVEVPGKVFRTEAFF;TRBV29-1;TRBJ1-1
clonotype100    CAARPGADKLIF;TRAV29/DV5;TRAJ34          CASSLEASGHPYEQYF;TRBV7-9;TRBJ2-7
clonotype101    CAVAIIIQGAQKLVF;TRAV41;TRAJ54           CASSLVGWKTSGFANTGELFF;TRBV11-1;TRBJ2-2
"""
# clone1 contains duplicates of the clonotypes that we wish to eliminate:
# clone1[clone1.index.duplicated(keep=False)] 
# OBS! Do something with these excluded clonotypes? Impute from these if first imputation fails?
"""
                genes_TRA                               genes_TRB
clonotype_id
clonotype10     CAASEPAKIPKAAGNKLTF;TRAV29/DV5;TRAJ17   CSVEVPGKVFRTEAFF;TRBV29-1;TRBJ1-1
clonotype10     CAVETSGSRLTF;TRAV39;TRAJ58              CSVEVPGKVFRTEAFF;TRBV29-1;TRBJ1-1
clonotype109    CAVGALGGTASKLTF;TRAV8-3;TRAJ44          CASRGIEQYF;TRBV27;TRBJ2-7
clonotype109    CVVLSGKLIF;TRAV24;TRAJ37                CASRGIEQYF;TRBV27;TRBJ2-7
"""
# remove clonotypes with multiple genes from imputation dict
clone2 = clone1.loc[clone1.index.drop_duplicates(keep=False)]
"""
                genes_TRA                       genes_TRB
clonotype_id
clonotype1      CAASQNEKLTF;TRAV21;TRAJ48       CASRIGAAGNSPLHF;TRBV27;TRBJ1-6
clonotype100    CAARPGADKLIF;TRAV29/DV5;TRAJ34  CASSLEASGHPYEQYF;TRBV7-9;TRBJ2-7
clonotype101    CAVAIIIQGAQKLVF;TRAV41;TRAJ54   CASSLVGWKTSGFANTGELFF;TRBV11-1;TRBJ2-2
clonotype102    CAVSDAGKSTF;TRAV1-2;TRAJ27      NaN
clonotype103    CAVGAGGTSYGKLTF;TRAV3;TRAJ52    CSGTITEQFF;TRBV29-1;TRBJ2-1
"""

#######################
## Clonotype duplicates
# List the clonotypes that share gene annotations but are named different
clone_dup = clone2[clone2.duplicated(subset=['genes_TRA','genes_TRB'], keep=False)].reset_index(drop=False).fillna('')

"""
    clonotype_id    genes_TRA                       genes_TRB
12  clonotype16     CIVRVGGDSWGKLQF;TRAV26-1;TRAJ24 CASSSLNTEAFF;TRBV27;TRBJ1-1
53  clonotype47     CAVRDTDARLMF;TRAV3;TRAJ31       CASSSVNEQYF;TRBV12-4;TRBJ2-7
0   clonotype104                                    CAWPRSTGELFF;TRBV30;TRBJ2-2
1   clonotype105                                    CAWPRSTGELFF;TRBV30;TRBJ2-2
2   clonotype116                                    CATRQNTEAFF;TRBV10-2;TRBJ1-1
3   clonotype117                                    CATRQNTEAFF;TRBV10-2;TRBJ1-1
4   clonotype125                                    CASRRAGPVPFF;TRBV19;TRBJ1-1
5   clonotype126                                    CASRRAGPVPFF;TRBV19;TRBJ1-1
6   clonotype130                                    CASSDLNSPLHF;TRBV27;TRBJ1-6
7   clonotype131                                    CASSDLNSPLHF;TRBV27;TRBJ1-6
8   clonotype141   CIVRVGGDSWGKLQF;TRAV26-1;TRAJ24  CASSSLNTEAFF;TRBV27;TRBJ1-1
9   clonotype142   CIVRVGGDSWGKLQF;TRAV26-1;TRAJ24  CASSSLNTEAFF;TRBV27;TRBJ1-1
10  clonotype143   CIVRVGGDSWGKLQF;TRAV26-1;TRAJ24  CASSSLNTEAFF;TRBV27;TRBJ1-1
"""

# Sort df by the numerical value of clonotype id: clonotyp438
clone_dup.sort_values('clonotype_id', key = lambda x: (x.str.split('clonotype').str[-1].astype(int)), inplace=True)

# 1. Dict of all duplicated clones. 
# Key is the clonotype id of a duplicated clone. Value is the first instance of that clone
clone_repl = dict()
def duplicated_clones(x, dct):
    for i in x[1:]:
        dct[i] = x[0]
        
"""
{'clonotype126': 'clonotype125',
 'clonotype393': 'clonotype392',
 'clonotype217': 'clonotype216',
 'clonotype218': 'clonotype216',
 'clonotype131': 'clonotype130',
 'clonotype379': 'clonotype378',
 'clonotype380': 'clonotype378',
 'clonotype396': 'clonotype395',
 'clonotype397': 'clonotype395',
 'clonotype165': 'clonotype164',
 'clonotype468': 'clonotype467',
 'clonotype470': 'clonotype469',
 'clonotype471': 'clonotype469',
 'clonotype390': 'clonotype389',
 'clonotype391': 'clonotype389',
 'clonotype223': 'clonotype222',
 'clonotype464': 'clonotype463',
 'clonotype465': 'clonotype463',
 'clonotype473': 'clonotype472',
 'clonotype474': 'clonotype472',
 'clonotype294': 'clonotype293',
 'clonotype401': 'clonotype400',
 'clonotype402': 'clonotype400',
 'clonotype227': 'clonotype226',
 'clonotype386': 'clonotype385',
 'clonotype205': 'clonotype204',
 'clonotype310': 'clonotype309',
 'clonotype311': 'clonotype309',
 'clonotype117': 'clonotype116',
 'clonotype105': 'clonotype104',
 'clonotype210': 'clonotype209',
 'clonotype146': 'clonotype47',
 'clonotype141': 'clonotype16',
 'clonotype142': 'clonotype16',
 'clonotype143': 'clonotype16'}
"""
# Populate clone_repl
clone_dup.groupby(['genes_TRA','genes_TRB']).clonotype_id.unique().apply(lambda x: duplicated_clones(x, clone_repl))
# Annotate the duplicated clonotypes with the first occuring clonotype in order to replace them.
df['clone_repl'] = df.clonotype.map(clone_repl)

####################################
# 2. Dict of unique clonotype labels

# drop clonotypes that have identical gene annotation (keep first)
clone2.drop_duplicates(subset=['genes_TRA','genes_TRB'], inplace=True)
# Dict keys: genes_TRA & genes_TRB, values: the first clonotype id
clone_dct = clone2.reset_index(drop=False).fillna('').groupby(['genes_TRA','genes_TRB']).clonotype_id.unique().apply(lambda x: x[0]).to_frame()


###########################################
#          Impute NaN clonotypes          #
###########################################
print('impute nan clonotypes')

df['imp_clone'] = np.nan
df['imp_type'] = np.nan

# For all the GEMs that did not have a clonotype annotation, search the clonotype dictionary
for gem, grp in df[df.clonotype.isnull()].groupby(['gem']):
    chain_a = grp[grp.chain == 'TRA']
    chain_b = grp[grp.chain == 'TRB']
    
    # List TRA annotations
    if len(chain_a) > 0:
        gene_a = [g for g in chain_a.genes]
    # If no TRA annotations exist: find the matching TRA's to the relevant TRBs
    else:
        gene_a = list()
        for g in chain_b.genes:
            dct_bool = clone_dct.index.isin([g], level=1)
            if dct_bool.any():
                # List all A-chains related to the given B-chain
                dct_idx = clone_dct[dct_bool].index.get_level_values(0)[0]
                gene_a.append(dct_idx)
                
    # List TRB annotations
    if len(chain_b) > 0:
        gene_b = [g for g in chain_b.genes]
    else:
        gene_b = list()
        for g in chain_a.genes:
            dct_bool = clone_dct.index.isin([g], level=0)
            if dct_bool.any():
                # List all B-chains related to the given A-chain
                dct_idx = clone_dct[dct_bool].index.get_level_values(1)[0]
                gene_b.append(dct_idx)
    
    # List the genes that are found in dictionary of clones
    genes = [g for g in itertools.product(gene_a, gene_b) if clone_dct.index.isin([g]).any()]
    
    # Extract the relevant entries from dict
    tmp_dct = clone_dct.loc[genes, 'clonotype_id']
    """
    genes_TRA                     genes_TRB                       
    CALRDMEYGNKLVF;TRAV16;TRAJ47  CASSLIVSGGANEQFF;TRBV7-3;TRBJ2-1    clonotype2
    CAVNYGNKLVF;TRAV3;TRAJ47      CASSLVAGAPSEQYF;TRBV27;TRBJ2-7      clonotype3
    """
    # Convert the multiindex to single index
    tmp_dct.index = tmp_dct.index.to_flat_index()
    """
    (CALRDMEYGNKLVF;TRAV16;TRAJ47, CASSLIVSGGANEQFF;TRBV7-3;TRBJ2-1)    clonotype2
    (CAVNYGNKLVF;TRAV3;TRAJ47, CASSLVAGAPSEQYF;TRBV27;TRBJ2-7)          clonotype3
    """
    # Convert df to contain a chain per row (hence 2 rows per clonotype id)
    tmp_dct = tmp_dct.reset_index(drop=False).explode('index').set_index('index')
    """
                                     clonotype_id
    index                                        
    CALRDMEYGNKLVF;TRAV16;TRAJ47       clonotype2
    CASSLIVSGGANEQFF;TRBV7-3;TRBJ2-1   clonotype2
    CAVNYGNKLVF;TRAV3;TRAJ47           clonotype3
    CASSLVAGAPSEQYF;TRBV27;TRBJ2-7     clonotype3
    """
    # Remove duplicated genes (if a geneset is annotated with multiple times with different clonotypes)
    tmp_dct = tmp_dct[~tmp_dct.index.duplicated(keep='first')]
    
    # Impute: Map clonotypes by genes onto DF
    df.loc[grp.index, 'imp_clone'] = df.loc[grp.index].set_index('genes').index.map(tmp_dct.clonotype_id)

    # Annotate type of imputation
    if tmp_dct.clonotype_id.nunique() > 1:
        df.loc[grp.index, 'imp_type'] = 'multiplet'
    elif tmp_dct.clonotype_id.nunique() == 1:
        df.loc[grp.index, 'imp_type'] = 'singlet'
        
    if df.loc[grp.index].imp_clone.isna().all():
        df.loc[grp.index, 'imp_type'] = 'new'
    elif df.loc[grp.index].imp_clone.isna().any():
        df.loc[grp.index, 'imp_type'] = 'multiplet'

##################################
## Name non-consensus clonotypes 
# 1. combine all genes of each GEM into a single variable
# 2. then group all GEMs by this variable to form new clonotypes
# 3. use the index to name the new clonotypes
new_clones_dct = 'new' + (df[df.imp_type == 'new'].groupby(['gem']).genes.unique().apply(lambda x: ' | '.join(x))
                          .to_frame().reset_index(drop=False).groupby('genes')['gem'].unique().to_frame()
                          .reset_index(drop=False).reset_index(drop=False).rename(columns={'index':'imp_clone'})
                          .explode('gem').set_index('gem').imp_clone.astype(str))
"""
gem
GACGCGTTCGTAGGTT-1      new0
TTCTCAATCGTCACGG-1      new0
AGTAGTCCAGACGCAA-1      new1
GAGCAGAGTAAGTTCC-1      new2
CGTTGGGTCGTCCGTT-1      new3
"""
df['new_clone'] = df.gem.map(new_clones_dct)

###########################################################
## Rank the chains of the clonotypes by UMI (10x and new) 
new_clones_df = df[df.imp_type == 'new'].sort_values(by=['gem','genes'])
# 1. Sum UMI per chain per clonotype
# 2. Sort the chains within a clonotype per UMI sum
# 3. Use this sorting as a ranking of chains within a clonotype
rank_chain_dct_new = (new_clones_df.groupby(['new_clone','genes']).umi_count.sum().to_frame()
                      .reset_index(drop=False).sort_values(by=['new_clone','umi_count'])
                      .reset_index(drop=True).reset_index(drop=False).set_index(['new_clone','genes'])['index'])

"""
new_clone  genes                            
new0       CAAENQGGKLIF;TRAV29/DV5;TRAJ23         0
new1       CAAFRAYNTDKLIF;TRAV29/DV5;TRAJ34       1
new10      CAASQNEKLTF;TRAV21;TRAJ48              2
           CASSLIVSGGANEQFF;TRBV7-3;TRBJ2-1       3
new100     CAVETSGSRLTF;TRAV39;TRAJ58             4
"""
df['chain_rank_new'] = df.set_index(['new_clone','genes']).index.map(rank_chain_dct_new)

################################################################
#                    Compile all clonotypes                    #
################################################################
max_ct_num = df.clonotype.dropna().str.split('clonotype').str[1].dropna().astype(int).max()
def annotate_clonotype(row, max_ct_num=max_ct_num):
    if row.clone_repl == row.clone_repl:
        return int(row.clone_repl.split('clonotype')[1])
    elif row.imp_type == 'singlet':
        return int(row.imp_clone.split('clonotype')[1])
    elif row.new_clone == row.new_clone:
        return int(row.new_clone.split('new')[1]) + max_ct_num + 1
    elif row.clonotype == row.clonotype:
        spl = row.clonotype.split('clonotype')
        if len(spl) == 1:
            return 0 #np.nan
        return int(row.clonotype.split('clonotype')[1])
    else:
        print('No-annotation clonotypes', row.clonotype, row.gem)

#################################################################
#        Filter - Remove GEMs with no credible clonotype        #
#################################################################
flt_df = df.dropna(subset=['clonotype','imp_clone','new_clone'], how='all').copy()
flt_df = flt_df[flt_df.imp_type != 'multiplet']

flt_df['ct'] = flt_df.apply(lambda row: annotate_clonotype(row), axis=1) # OBS! changed from ct to ct_sub
assert flt_df.ct.isna().any() == False
################################################################
#                       Augment by chain                       #
################################################################
flt_df['chain_count'] = flt_df.groupby(['gem', 'chain']).contig_id.transform('size')

tra_df = flt_df[(flt_df.chain == "TRA")].copy()
trb_df = flt_df[(flt_df.chain == "TRB")].copy()

rprt['tra_gems'] = list(set(tra_df.gem))
rprt['trb_gems'] = list(set(trb_df.gem))

# Group all GEMs, sort chains by increasing UMI count and if no difference then sort by accumulated umi per clonotype
tra_df.sort_values(by=['gem', 'umi_count', 'is_cell', 'high_confidence'], inplace=True)
trb_df.sort_values(by=['gem', 'umi_count', 'is_cell', 'high_confidence'], inplace=True)

# Replacing genes annotated with None with ''
tra_df['genes'] = tra_df.replace([None], ['']).apply(lambda x: ';'.join(x[['v_gene','j_gene','c_gene']]), axis=1)
trb_df['genes'] = trb_df.replace([None], ['']).apply(lambda x: ';'.join(x[['v_gene','d_gene','j_gene','c_gene']]), axis=1)

def annotate_umi_lst(df):
    dct = df.groupby(['gem','cdr3']).umi_count.apply(sum).to_frame().groupby('gem').umi_count.apply(lambda x: sorted(np.array(x))).to_dict()
    return df.gem.map(dct)

def annotate_cdr3_lst(df):
    dct = df.groupby(['gem']).cdr3.unique().to_dict()
    return df.gem.map(dct)

def annotate_genes_lst(df):
    # I use unique instead of listing because the same CDR3 can occur multiple times in a GEM.
    dct = df.groupby(['gem']).genes.unique().to_dict()
    return df.gem.map(dct)

def annotate_single(df):
    return df.umi_count_lst.apply(lambda x: True if len(x)==1 else False)

def annotate_delta_umi(df):
    def calc_delta(x):
        if len(x) == 1:
            return x[-1]/0.25
        elif len(x) == 0:
            return 0
        else:
            return (x[-1])/(x[-2]+0.25)
    return df.umi_count_lst.apply(calc_delta)

tra_df['umi_count_lst'] = annotate_umi_lst(tra_df)
trb_df['umi_count_lst'] = annotate_umi_lst(trb_df)

tra_df['cdr3_lst'] = annotate_cdr3_lst(tra_df)
trb_df['cdr3_lst'] = annotate_cdr3_lst(trb_df)

tra_df['genes_lst'] = annotate_genes_lst(tra_df)
trb_df['genes_lst'] = annotate_genes_lst(trb_df)

tra_df['single'] = annotate_single(tra_df)
trb_df['single'] = annotate_single(trb_df)

tra_df['delta_umi'] = annotate_delta_umi(tra_df)
trb_df['delta_umi'] = annotate_delta_umi(trb_df)

# Assert that the annotated CDR3 lst is as long as the UMI lst
assert tra_df.apply(lambda x: True if len(x.umi_count_lst) == len(x.cdr3_lst) else False, axis=1).all()

# ### Keep chain with highest UMI count and merge
tra_df.drop_duplicates(subset=['gem'], keep='last', inplace=True)
trb_df.drop_duplicates(subset=['gem'], keep='last', inplace=True)

###############################################################
#                            Merge                            #
###############################################################
# In case merge doesn't accept NaNs in the key of merging
tra_df.clonotype.fillna('', inplace=True)
trb_df.clonotype.fillna('', inplace=True)

# It is not necessary to merge on clonotype, only to keep a single column of clonotype instead of two
tcr_df = pd.merge(tra_df, trb_df, how='outer', on=['gem','ct','clonotype'], suffixes=('_TRA', '_TRB'))
tcr_df.clonotype.replace('',np.nan, inplace=True)

assert tcr_df.shape[0] == tcr_df.gem.unique().shape[0]

################################################################
#                 Define functional clonotypes                 #
################################################################
# Functional clonotypes are based on only CDR3
# Get (if available) both CDR3 sequences defined per sub clonotype
# Sort the CDR3 pairs by frequency, to name highly frequent CDR3 pairs with a low functional clonotype ID
# Is not working and therefore not implemented...
func = tcr_df.groupby(['ct','cdr3_TRA','cdr3_TRB']).size().to_frame().reset_index().drop_duplicates(subset='ct', keep='last').rename(columns={0:'ct_freq'})
func.sort_values(by='ct_freq', ascending=False, inplace=True)
func.reset_index(drop=True, inplace=True)
func.reset_index(inplace=True)
func.rename(columns={'index':'ct_super'}, inplace=True)
func.ct_super = func.ct_super + 1
dct = func.set_index('ct').ct_super
tcr_df['ct_super'] = tcr_df.set_index('ct').index.map(dct)

################################################################
#                        Augment by GEM                        #
################################################################
tcr_df['num_clonotype'] = pd.to_numeric(tcr_df['clonotype'].fillna('None').str.split('clonotype').str[1],
                                        errors='coerce').replace(np.nan, 0, regex=True).astype(int)

tcr_df['single_chain_only'] = tcr_df[['chain_TRA', 'chain_TRB']].isna().any(axis=1)
tcr_df['umi_count_tcr'] = tcr_df.umi_count_TRA.fillna(0) + tcr_df.umi_count_TRB.fillna(0)
tcr_df['cdr3_comb'] = tcr_df.cdr3_TRA.fillna('') + tcr_df.cdr3_TRB.fillna('')

def define_tcr_categories(row):
    if (row['single_TRA'] == True) & (row['single_TRB'] == True):
        return 'unique chains'
    if row['single_chain_only']:
        return 'missing chain'
    else:
        return 'multiple chains'

tcr_df['tcr_category'] = tcr_df.apply(lambda row: define_tcr_categories(row), axis=1)
tcr_df['no_filtration'] = True
tcr_df['exclude_single-chain_TCRs'] = tcr_df.apply(lambda row: True if (row.tcr_category == 'unique chains') or (row.tcr_category == 'multiple chains') else False, axis=1)
tcr_df['exclude_ambiguous_and_single-chain_TCRs'] = np.where((tcr_df.tcr_category == 'unique chains'), True, False)
tcr_df['exclude_ambiguous_TCRs'] = tcr_df.apply(lambda row: True if (row.tcr_category == 'unique chains') or (row.tcr_category == 'missing chains') else False, axis=1)
tcr_df['cell_flag'] = np.where(tcr_df.is_cell_TRA | tcr_df.is_cell_TRB, True, False)
tcr_df['cell_high_confidence'] = np.where(tcr_df.high_confidence_TRA & tcr_df.high_confidence_TRB, True, False)

# Label GEMs by which sorting (positive/negative)
if ('label_TRA' in tcr_df.columns):
    tcr_df['label'] = tcr_df.label_TRA
else:
    # Assuming all GEMs are from a positive sorting
    tcr_df['label'] = 1

assert tcr_df.ct.isna().any() == False


# ## Write data
new_column_order = ['gem', 'clonotype', 'num_clonotype', 'ct', #'ct_super',
                    'genes_TRA', 'genes_TRB', 'genes_lst_TRA', 'genes_lst_TRB',
                    'length_TRA', #'cdr1_TRA', 'cdr2_TRA',
                    'cdr3_TRA', 'umi_count_TRA', 'umi_count_lst_TRA', 'delta_umi_TRA', 'cdr3_lst_TRA',
                    'chain_count_TRA','single_TRA',
                    'length_TRB', #'cdr1_TRB', 'cdr2_TRB',
                    'cdr3_TRB', 'umi_count_TRB', 'umi_count_lst_TRB', 'delta_umi_TRB', 'cdr3_lst_TRB',
                    'chain_count_TRB','single_TRB',
                    'single_chain_only', 'umi_count_tcr', 'cdr3_comb', 'v_gene_TRA','j_gene_TRA', 'v_gene_TRB','j_gene_TRB',
                    'tcr_category', 'no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs','exclude_ambiguous_TCRs',
                    'label','cell_flag','cell_high_confidence']


tcr_df[new_column_order].to_csv(OUTPUT, index=False)

if report:
    with open(report, 'w') as outfile:
        json.dump(rprt, outfile)


