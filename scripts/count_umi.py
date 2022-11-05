#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np

MAP = snakemake.input.mpp #'experiments/exp13/run1/mhc/mapping/split/map/A1072B316' #
UMI = snakemake.input.umi #'experiments/exp13/run1/mhc/mapping/split/umi/A1072B316.tsv' #
OUTPUT = snakemake.output[0]

# Load
map_df = pd.read_csv(MAP, sep='\t', header=None, names=['seq','hits','barcode','read_id','label'], usecols=[0,1,5,6,8])
umi_df = pd.read_csv(UMI, sep='\t', header=None, names=['barcode','read_ext','umi'], skiprows=1, dtype=str)

map_df['read_id'], map_df['gem'] = map_df.read_id.str.split(r'\s+BX:Z:').str
map_df['read_ext'] = map_df.read_id.str.rsplit(':', n=1, expand=True)[1]
umi_df['read_ext'] = umi_df.read_ext.astype(str)
map_df['barcode'] = map_df.barcode.str.split(' ', expand=True)[0]

# To make sure that the "best" read is kept one could sort on hits, so we always keep the read with lowest hits.
map_df.drop_duplicates(subset=['read_id','gem'], keep='first', inplace=True)
map_df.drop_duplicates(subset=['read_id'], keep=False, inplace=True)

# To make sure that the "best" read is kept one could sort on umi length, so we always keep the read with full umi seq.
umi_df.replace('', np.nan, inplace=True)
umi_df.dropna(inplace=True)
umi_df.drop_duplicates(subset=['read_ext','barcode'], keep='first', inplace=True)
umi_df.drop_duplicates(subset=['read_ext'], keep=False, inplace=True)

df = pd.merge(map_df, umi_df, on=['read_ext','barcode'])

assert df.duplicated().sum() == 0

cnt = df.groupby(['gem','barcode']).umi.nunique().to_frame()
cnt['label'] = df.groupby(['gem','barcode']).label.unique().apply(lambda x: x[0] if len(x)==1 else np.nan)
cnt.to_csv(OUTPUT, header=False)

