#!/usr/bin/env python

import numpy as np
import pandas as pd
import re

# ######################
# #    SCRIPT THREE    #
# ######################
import numpy as np
import pandas as pd

# Combine chunks and chains

INPUTS = snakemake.input
OUTPUT = snakemake.output[0]

dfs = list()
for filename in INPUTS:
    tmp = pd.read_csv(filename, index_col=0)
    dfs.append(tmp)
    
df = pd.concat(dfs, axis=1)
df.loc['missing'] = 0 #Similarity score is 0 when CDR3 seq is ''
df['missing'] = 0

assert sum(df.index.duplicated()) == 0
assert sum(df.columns.duplicated()) == 0

df.to_csv(OUTPUT, index=True)

print(df)