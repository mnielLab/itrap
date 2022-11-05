#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')

#FILE = "/Volumes/tuba/herpov/tcr-pmhc-sc-project/data/exp5_MHC_IONTORRENT_archived/reports/read_counts_mhc_pipeline.txt"

df = pd.read_csv(snakemake.input[0], sep='\t')

df.rename(columns={'count':'frequency'}, inplace=True)
df["file"] = df.file.str.replace("cut_adapter", "Adapter trimming")
df["file"] = df.file.str.replace("outs", "GEM annotation")
df["file"] = df.file.str.replace("longranger_clean", "GEM filtering")
df["file"] = df.file.str.replace("output", "KMA mapping")

plt.bar(df.file.values, df.frequency.values)
plt.title("Read counts throughout pipeline")
plt.ylabel("Counts")
plt.xlabel("Pipeline steps")
plt.xticks(rotation=45, ha='right')
plt.savefig(snakemake.output[0], bbox_inches='tight')

