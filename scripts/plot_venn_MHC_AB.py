#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn2

# # Load
df = pd.read_csv(snakemake.input[0])

gem_mhc = df[df.umi_count_mhc.notnull()].gem
gem_cd8 = df[df.umi_count_cd8.notnull()].gem


# # Plot
a, b = set(gem_mhc), set(gem_cd8)

ab = len(a.intersection(b))
a_b = len(a) - ab
b_a = len(b) - ab

title = "GEMs containing MHC and/or Ab BCs"
venn2(subsets = (a_b, b_a, ab), set_labels = ('MHC', 'Ab'))
plt.title(title)
plt.savefig(snakemake.output[0], bbox_inches='tight')
