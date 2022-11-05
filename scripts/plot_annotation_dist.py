#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('ggplot')


# # Input
TCR_CL = snakemake.input[2] #'/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp9_TCR/augmented/tcr.clean.augmented.csv'
BC_MP = snakemake.input[0] #'/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp9_MHC_IONTORRENT/mapping/KMA-1t1/output/mapping.frag.gz'
BC_CL = snakemake.input[1] #'/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp9_MHC_IONTORRENT/mapping/KMA-1t1/output/mapping.clean.AKB.augmented.gz'
SPECIFICITY = snakemake.input[3] #'/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp9_MHC_IONTORRENT/barcode_library/barcode_specificity_annotations.xlsx'

# # Output

# # Load
sp_df = pd.read_excel(SPECIFICITY, sheet_name='CDX')
tcr_cl = pd.read_csv(TCR_CL)
bc_mp = pd.read_csv(BC_MP, sep="\t", usecols=[5,6], names=['barcode','read_id']) #, engine='python'
bc_mp['barcode'] = bc_mp.barcode.str.split(' ', expand=True)[0]
bc_mp['gem'] = bc_mp.read_id.str.split('BX:Z:', expand=True)[1]

bc_cl = pd.read_csv(BC_CL)


# # Set variables
a = set(tcr_cl.gem)
b = set(bc_cl.gem)
ab_bc = sp_df.Barcode

bc_cl['barcode_mhc'] = bc_cl.template_id_mhc.str.split('_', expand=True)[0]
bc_cl['barcode_cdx'] = bc_cl.template_id_cd8.str.split('_', expand=True)[0]


# # Plotting GEMs per barcode
data_structure = [('Mapped GEMs', (bc_mp[bc_mp.barcode.isin(ab_bc)].groupby('barcode').gem.unique().apply(len),
                                        bc_mp[~bc_mp.barcode.isin(ab_bc)].groupby('barcode').gem.unique().apply(len))),
                  ('Credibly mapped GEMs', (bc_cl.groupby('barcode_cdx').gem.unique().apply(len),
                                                 bc_cl.groupby('barcode_mhc').gem.unique().apply(len))),
                  ('Credibly mapped GEMs intersecting with TCR annotated GEMs', (bc_cl[bc_cl.gem.isin(a.intersection(b))].groupby('barcode_cdx').gem.unique().apply(len),
                                                                                      bc_cl[bc_cl.gem.isin(a.intersection(b))].groupby('barcode_mhc').gem.unique().apply(len))),
                  ('Credibly mapped GEMs not intersecting with TCR annotated GEMs', (bc_cl[bc_cl.gem.isin(b-a)].groupby('barcode_cdx').gem.unique().apply(len),
                                                                                          bc_cl[bc_cl.gem.isin(b-a)].groupby('barcode_mhc').gem.unique().apply(len)))]


fig, axes = plt.subplots(4,1, figsize=(10,14), sharex=True)

for ax, (title, dataset) in zip(axes, data_structure):
    for label, data in zip(['cdx','mhc'], dataset):
        x = data.index
        y = data.values
        ax.bar(x, y, label=label + ' (%d GEMs)'%sum(y))
        
    ax.set_xlim(-0.5, len(x)-0.5)
    ax.set_ylabel("GEM counts")
    ax.set_title(title)
    ax.legend()
plt.xlabel("Barcode templates")
plt.xticks(rotation=90, fontsize=7)
plt.savefig(snakemake.output[0], bbox_inches='tight')


# # Plotting GEMs per clonotype
tcr_data = [('GEMs with annotated clonotypes (%d GEMs)', tcr_cl.groupby('ct').gem.unique().apply(len)),
            ('Gems with annotated CT intersecting with\nGEMs having credibly mapped BCs (%d GEMs)', tcr_cl[tcr_cl.gem.isin(a.intersection(b))].groupby('ct').gem.unique().apply(len)),
            ('Gems with annotated CT not intersecting with\nGEMs having credibly mapped BCs (%d GEMs)', tcr_cl[tcr_cl.gem.isin(a - b)].groupby('ct').gem.unique().apply(len))]

fig, axes = plt.subplots(3,1, figsize=(10,14), sharex=True)

for ax, (title, data) in zip(axes, tcr_data):
    x = data.index
    y = data.values
    ax.bar(x, y)
        
    ax.set_xlim(-0.5, len(x)-0.5)
    ax.set_ylabel("GEM counts")
    ax.set_title(title %sum(y))

plt.xlabel("Clonotypes in numerical order")
plt.xticks(rotation=90, fontsize=7)
plt.savefig(snakemake.output[1], bbox_inches='tight')

