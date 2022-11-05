#!/opt/anaconda3/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import math

import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects

plt.style.use('ggplot')

# Input
#KMA-1t1/output/mapping.clean.gz

# Load data
df = pd.read_csv(snakemake.input[0]) 

# Main

# Uncertainty distribution
plt.bar(df.uncertainty.value_counts().index, df.uncertainty.value_counts().values)
plt.title("Distribution of uncertainties in mapping")
plt.xlabel('Uncertainty')
plt.ylabel('Read counts')
plt.savefig(snakemake.output[0], bbox_inches='tight') #OUT_DIR + "uncertainty_dist.pdf"

# Distribution of reads per template
x = df.groupby('barcode').size().index
y = df.groupby('barcode').size().values
plt.figure(figsize=(20,7))
plt.bar(x,y)
plt.xlim(-0.5, len(x)-0.5)
plt.xticks(rotation=90, fontsize=8)
plt.xlabel("Barcode templates")
plt.ylabel("Read counts")
plt.title("Read count per barcode template")
plt.savefig(snakemake.output[1], bbox_inches='tight') #OUT_DIR + "reads_per_template.pdf"

# Distribution of GEMs per template
x = df.groupby(['barcode']).gem.unique().apply(len).index
y = df.groupby(['barcode']).gem.unique().apply(len).values
plt.figure(figsize=(20,7))
plt.bar(x,y)
plt.xlim(-0.5, len(x)-0.5)
plt.xticks(rotation=90, fontsize=8)
plt.xlabel("Barcode templates")
plt.ylabel("GEM counts")
plt.title("GEM count per barcode template")
plt.savefig(snakemake.output[2], bbox_inches='tight') #OUT_DIR + "gems_per_template.pdf"

# Reads per GEM
x = df.groupby(['gem']).query_id.unique().apply(len).index
y = df.groupby(['gem']).query_id.unique().apply(len).values
plt.figure(figsize=(20,7))
plt.bar(x,y)
plt.xlim(-0.5, len(x)-0.5)
plt.tick_params(labelbottom=False)
plt.xlabel("GEMs")
plt.ylabel("Read counts")
plt.title("Read count per GEM")
plt.savefig(snakemake.output[3], bbox_inches='tight') #OUT_DIR + "reads_per_gem.pdf"

# Read alignments
cc = plt.rcParams['axes.prop_cycle'].by_key()['color'] * math.ceil(df.uncertainty.max()/7)
data_categories = {'expected_templates':{'TSO':(0, 13, cc[3]),
                                         'Primer_B':(13, 19, cc[0]),
                                         'N6':(32, 6, cc[4]),
                                         'Oligo_B':(38, 25, cc[2]),
                                         'Anneal':(63, 25, cc[1]),
                                         'Oligo_A':(88, 25, cc[2]),
                                         'N6-':(113, 6, cc[4]),
                                         'Primer_A':(119, 18, cc[0]),
                                         'ID':(137, 8, cc[5])}}

for data_category in data_categories:
    title = "Alignment summary"
    # Set up the axes with gridspec
    fig = plt.figure(figsize=(14, 6))
    grid = plt.GridSpec(2, 1, hspace=0.25, height_ratios=[10,1])
    # Main
    main_ax = plt.subplot(grid[0, 0])
    
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.viridis(np.linspace(0,1,len(df.uncertainty.unique()))))
    cc2 = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i, u in enumerate(sorted(df.uncertainty.unique())):
        main = plt.plot([df[df.uncertainty == u].t_alignment_start.to_list(), df[df.uncertainty == u].t_alignment_end.to_list()],
                        [df[df.uncertainty == u].score.to_list(),df[df.uncertainty == u].score.to_list()], '.-', linewidth=0.2, c=cc2[i])
    main = plt.xlim(-3,149) # Variable!
    main = plt.xlabel("Alignment positions in template")
    main = plt.ylabel("Alignment score")
    main = plt.title(title)
    # Add legend
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=1, vmax=df.uncertainty.max()))
    if mpl.__version__ < '3.1':
        sm._A = []
    cbaxes = fig.add_axes([0.9, 0.27, 0.01, 0.61]) #[x-pos,y-pos,width, height]
    cbar = plt.colorbar(sm, cax = cbaxes) #pad=-0.0005
    cbar.set_ticks(df.uncertainty.unique())
    cbar.set_ticklabels(df.uncertainty.unique())
    cbar.ax.set_ylabel('# of equally likely template matches', rotation=90) #270
    # Subs
    subs_ax = plt.subplot(grid[1, 0], yticklabels=[], yticks=[], sharex=main_ax)
    # Add template design
    for k, v in data_categories[data_category].items():
        print(k,v)
        # Create a Rectangle patch
        rect = patches.Rectangle((v[0],0),v[1],1, fill=True, color=v[2]) #linewidth=1,edgecolor='r',facecolor='none'
        subs_ax.add_patch(rect)
        # Add text
        block = " ".join(k.split("_")).split("-")[0]
        txt = subs_ax.text(v[1]/2.0+v[0], 0.5, block, horizontalalignment='center', verticalalignment='center', color='white',
                           size=15, weight='heavy')
        txt.set_path_effects([PathEffects.withStroke(linewidth=0.5, foreground='black')])       
    plt.savefig(snakemake.output[4], bbox_inches='tight') #OUT_DIR + "alignments.png"

