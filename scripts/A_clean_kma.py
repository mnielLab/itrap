#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd


kma_fil = snakemake.input[0] #"/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC_" + PLATFORM + "/mapping/kma-1t1/output/mapping.frag.gz"
clean_fil = snakemake.output[0] #"/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC_" + PLATFORM + "/mapping/kma-1t1/output/mapping.clean.gz"


# # Main
df = pd.read_csv(kma_fil, sep='\t', usecols = (1, 2, 3, 4, 5, 6), names = ('uncertainty', 'score', 't_alignment_start', 't_alignment_end', 'template_id', 'read_header'))
df['query_id'], df['gem'] = df['read_header'].str.split(" ", n = 1).str
df.replace({"gem": r"^BX:Z:"}, {"gem": ""}, regex=True, inplace=True)
assert df.duplicated(['query_id']).sum() == 0

df['credible_alignment'] = np.where(df['uncertainty'] == 1, True, False)
df['alignment_length'] = df['t_alignment_end'] - df['t_alignment_start']
df['template_id'] = df.template_id.str.split(" ", n=1,expand=True)[0] + "_sample" # HACK!
df['barcode'], df['sample'] = df.template_id.str.rsplit("_", n=1).str

df.to_csv(clean_fil, index=False)


