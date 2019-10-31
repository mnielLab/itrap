#!/opt/anaconda3/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects

plt.style.use('ggplot')
cc = plt.rcParams['axes.prop_cycle'].by_key()['color'] * 2

# ARGS
EXP = "exp2"
PRJ = "blast"
PLATFORM = "_ILLUMINA"

# Plotting directory
FIG_DIR = "/home/tuba/herpov/tcr-pmhc-sc-project/plots/" + EXP + "/" + PRJ + "/"

# Input data
BLAST_DIR = "/home/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC" + PLATFORM + "/" + PRJ + "/" # OBS!
HEAD_FILE = "/home/tuba/herpov/tcr-pmhc-sc-project/data/" + EXP + "_MHC" + PLATFORM + "/cutadapt_naming/no_adapters/headers/annotated_headers.lst"

# Load header data
head_df = pd.read_csv(HEAD_FILE, sep=" ", names=['query_id', 'gem', 'tso', 'b_primer', 'anneal', 'cd8_primer', 'mhc_primer'])
head_df['a_primer'] = np.where(head_df.cd8_primer == 'no_adapter', head_df.mhc_primer, head_df.cd8_primer)
head_df.replace('no_adapter', np.NaN, inplace=True)
head_df['match'] = head_df.loc[:, ['tso', 'b_primer', 'anneal', 'a_primer']].count(axis=1)

# Main
data_categories = {'expected_templates':{'TSO':(0, 13, cc[3]), 'Primer_B':(13, 19, cc[0]),'N6':(32, 6, cc[4]),'Oligo_B':(38, 25, cc[2]),'Anneal':(63, 25, cc[1]),'Oligo_A':(88, 25, cc[2]),'N6-':(113, 6, cc[4]),'Primer_A':(119, 18, cc[0]),'ID':(137, 8, cc[5])},
                   'overrep_seq_templates':{'TSO':(0, 13, cc[3]), 'Some_sequence':(13, 70, cc[4])},
                   'reversed_templates':{'TSO':(0, 13, cc[3]),'ID':(13, 8, cc[5]),'Primer_A':(21, 18, cc[0]),'N6-':(39, 6, cc[4]),'Oligo_A':(45, 25, cc[2]),'Anneal':(70, 25, cc[1]),'Oligo_B':(95, 25, cc[2]),'N6':(120, 6, cc[4]),'Primer_B':(126, 19, cc[0])},
                   'rev_templates_start':{'TSO':(0, 13, cc[3]),'ID':(13, 8, cc[5]),'Primer_A':(21, 18, cc[0])}}



blast_hit_count, read_hits = dict(), dict()

for data_category in data_categories:
    print(data_category)
    
    BLAST_TABLE = BLAST_DIR + data_category + "/ILLUMINA-MHC-2.blast.tsv" # OBS! Changes!
    OUT_TABLE = BLAST_DIR + data_category + "/blast.annotated.tsv"
    
    # Load data
    blast_df = pd.read_csv(BLAST_TABLE, sep="\t", names=['query_id', 'template_id', 'identity', 'alignment_length',
                                                         'mismatches', 'gap_openings', 'q_alignment_start',
                                                         'q_alignment_end', 't_alignment_start', 't_alignment_end',
                                                         'e-value', 'bit_score'])
    
    # Get best BLAST hit
    uniq_blast_df = blast_df.drop_duplicates(['query_id', 'gem', 'bit_score', 't_alignment_start', 't_alignment_end', 'q_alignment_start', 'q_alignment_end', 'match'], keep=False).sort_values('bit_score', ascending=False).drop_duplicates(['query_id'])
    
    blast_hit_count[data_category] = uniq_blast_df.shape[0]
    
    read_hits[data_category] = set(uniq_blast_df.query_id.to_list())

    ann_df = pd.merge(blast_df, head_df[['query_id','gem','tso','b_primer','anneal','a_primer','match']], how='left', on=['query_id'])
    ann_df.to_csv(OUT_TABLE, sep="\t", index=False)
    
    
    # Sort BLAST hits according to unique templates
    template_counts_df = uniq_blast_df.set_index(['template_id', 'query_id']).count(level='template_id').sort_values('identity', ascending=False)
    
    ################################################################################################################
    ###                                                PLOTTING                                                  ###
    ################################################################################################################
    
    # Plot
    title = "Read counts for each BLAST template"
    fig = plt.figure(figsize=(14,6))
    plt.bar(x=template_counts_df.index, height=template_counts_df.identity)
    plt.title(title)
    plt.xlabel("Template IDs")
    plt.ylabel("Read counts")
    plt.xticks(rotation=90, fontsize=4) # ha='right', 
    plt.xlim(-0.5,len(template_counts_df.index)-.5) # Remove empty margins within plot
    
    plt.savefig(FIG_DIR + "_".join(title.split(" ")) + "_" + data_category + ".pdf", bbox_inches='tight')
    #plt.show()

    # Plot
    title = "Distributions of alignment metrics"
    fig, axs = plt.subplots(1, 4, figsize=(18, 6))
    axs[0].hist(uniq_blast_df.identity, bins=range(100+2))
    axs[1].hist(uniq_blast_df.alignment_length, bins=range(max(uniq_blast_df.alignment_length)+2))
    axs[2].hist(uniq_blast_df.mismatches, bins=range(max(uniq_blast_df.mismatches)+2))
    axs[3].hist(uniq_blast_df.gap_openings, bins=range(max(uniq_blast_df.gap_openings)+2))
    
    axs[0].set_ylabel("Frequency")
    axs[0].set_xlabel("Identity")
    axs[1].set_xlabel("Alignment length")
    axs[2].set_xlabel("Mismatches")
    axs[3].set_xlabel("Gap openings")
    fig.suptitle(title, fontsize=16)
    
    plt.savefig(FIG_DIR + "_".join(title.split(" ")) + "_" + data_category + ".pdf", bbox_inches='tight')
    #plt.show()

    # Plot
    title = "Alignment positions"
    # Set up the axes with gridspec
    fig = plt.figure(figsize=(14, 6))
    grid = plt.GridSpec(2, 1, hspace=0.25, height_ratios=[10,1])
    # Main
    main_ax = plt.subplot(grid[0, 0])

    for i in range(ann_df.match.max() + 1):
        main = plt.plot([ann_df[ann_df.match == i].t_alignment_start.to_list(), ann_df[ann_df.match == i].t_alignment_end.to_list()],
                        [ann_df[ann_df.match == i].q_alignment_start.to_list(),ann_df[ann_df.match == i].q_alignment_end.to_list()],
                        '.-', c=cc[i], label="%i block matches" %i)
    main = plt.xlim(-3,153) # Variable!
    main = plt.xlabel("Alignment positions in template")
    main = plt.ylabel("Alignment positions in query")
    main = plt.title(title)
    # Subs
    subs_ax = plt.subplot(grid[1, 0], yticklabels=[], sharex=main_ax)
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

    plt.savefig(FIG_DIR + "_".join(title.split(" ")) + "_" + data_category + ".png", bbox_inches='tight')
    #plt.show()

# Some plots have not been included in this script
