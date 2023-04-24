
"""
Run pipeline from location of the Snakefile:
snakemake --config exp=exp13 run=run1 --cores 4 --use-conda

Pipeline requires:
Input:
    - Cellranger multi output of TCR, and pMHC as a minimum.
    - Additional Cellranger features such as gene expression (GEX) and cell hashing barcodes can be handled.
    - An experimental design template (see barcode_specificity_annotations.xlsx).
"""
import os

EXPERIMENTAL_DESIGN_TEMPLATE = "experiments/exp13/run1/lib/barcode_specificity_annotations.xlsx"
CELLRANGER_DIR = "experiments/exp13/run1/tcr/cellranger_pos"

WRK_DIR = workflow.basedir
EXP_DIR = os.path.join(WRK_DIR, "experiments", config["exp"], config["run"])
RES_DIR = os.path.join(EXP_DIR, "res")
PLT_DIR = os.path.join(EXP_DIR, "plt")


rule all:
    input: EXP_DIR + '/done.ok'

#################################################################
#                         TCR clonotypes                        #
#################################################################
rule clean_augment_tcr:
    """
    Reannotating clonotypes:
    - Clonotypes of identical TCRab AA sequences are merged
    - GEMs with no clonotypes may be assigned an clonotype in two ways:
      1. if the TCRab pairs uniquely match an existing clonotype 
      2. if the TCRab are represent a novel clonotype
    """
    input:
        contig = CELLRANGER_DIR + "/outs/multi/vdj_t/all_contig_annotations.csv"
    params:
        clonot = CELLRANGER_DIR + "/outs/per_sample_outs/cellranger_pos/vdj_t/consensus_annotations.csv"
    output:
        output = RES_DIR + "/tables/tcr.clean.augmented.csv"
    conda:
        "envs/basic_dependencies.yaml"
    shell:
        "python scripts/A_clean_augment_tcr.py \
            --contig {input.contig} \
            --consensus {params.clonot} \
            --output {output.output}"
        

#################################################################
#               Reformat Cellranger barcode output              #
#################################################################
rule parse_count_matrix:
    """
    Run this rule if barcode reads origin from 10x barcodes and are mapped using Cellranger.
    Annotate types of barcode: pMHC, cell hashing, or surface markers.
    Annotate most abundant barcode per GEM and other summary metrics.
    """
    input:
        brc = CELLRANGER_DIR + "/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz",
        gem = CELLRANGER_DIR + "/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz",
        mtx = CELLRANGER_DIR + "/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz",
        lbl = EXPERIMENTAL_DESIGN_TEMPLATE,
        ann = WRK_DIR + "/tools/detected_responses_annotation.xlsx"
    output:
        RES_DIR + "/tables/brc.augmented.csv"
    conda:
        "envs/basic_dependencies.yaml"
    shell:
        "python scripts/A_parse_count_matrix.py \
            --features {input.brc} \
            --barcodes {input.gem} \
            --matrix {input.mtx} \
            --multimers {input.lbl} \
            --responses {input.ann} \
            --output {output}"
        

#################################################################
#                        Gene expression                        #
#################################################################
ruleorder: filter_gex_data > no_gex_data

rule filter_gex_data:
    """
    Filtering data based on Gene expression filters, using Seurat:
    https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
    """
    input:
        CELLRANGER_DIR + "/outs/multi/count/raw_feature_bc_matrix.h5"
    output:
        RES_DIR + "/tables/gex.txt",
        PLT_DIR + "/GEX/filtering/violin.png",
        PLT_DIR + "/GEX/filtering/scatter.png"
    conda:
        "envs/seurat_gex_filtering.yaml"
    shell:
        "Rscript ./scripts/B0_seurat_gex_analysis.R {input} {output}"
        
rule no_gex_data:
    """
    Generate a dummy GEX filtering file.
    """
    output:
        touch(RES_DIR + "/tables/gex.txt")
#################################################################
#             Seurat: Demultiplexing Hashtag Oligos             #
#################################################################
rule prep_seurat_HTO_analysis:
    """
    Convert count matrix to RDS format for Seurat analysis.
    """
    input:
        df = rules.clean_augment_tcr.output.output,
        hto = EXPERIMENTAL_DESIGN_TEMPLATE,
        brc = CELLRANGER_DIR + "/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz",
        gem = CELLRANGER_DIR + "/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz",
        mtx = CELLRANGER_DIR + "/outs/multi/count/raw_feature_bc_matrix/matrix.tsv.gz",
    output:
        rds = temp(RES_DIR + "/tables/hto_count_matrix.rds")
    conda:
        "envs/prep_seurat.yaml"
    script:
        "scripts/B1_seurat_hto_prep.py"

ruleorder: seurat_HTO_analysis > no_seurat_HTO_analysis

rule seurat_HTO_analysis:
    """
    Analyze cell hashing barcodes with Seurat HTO method.
    """
    input:
        rules.prep_seurat_HTO_analysis.output.rds
    output:
        out_file = RES_DIR + "/tables/hto.csv",
        ridge_plot = PLT_DIR + "/HTO/ridgeplot.png",
        violin_plot = PLT_DIR + "/HTO/violinplot.png",
        heatmap = PLT_DIR + "/HTO/heatmap.png"
    conda:
        "envs/seurat_hto.yaml"
    shell:
        "Rscript ./scripts/B2_seurat_hto_analysis.R \
            {input} \
            {output.out_file} \
            {output.ridge_plot} \
            {output.violin_plot} \
            {output.heatmap}"
 
rule no_seurat_HTO_analysis:
    """
    If no cell hashing barcodes are available, mimic the HTO output.
    """
    output:
        out_file = RES_DIR + "/tables/hto.csv"
    run:
        import pandas as pd
        
        df = pd.DataFrame(columns=['gem','seurat','umi_count_hto','feature_hto','hto_max_id','hto_sec_id',
                                   'hto_margin','hto_classification','hto_global_class','hash_id'])
        df.to_csv(output.out_file, index=False)
        

#################################################################
#                       Concatenate data                        #
#################################################################
rule comb_barcodes_TCR:
    """
    Merge barcode and TCR data based on GEM barcode.
    """
    input:
        tcr = rules.clean_augment_tcr.output.output,
        brc = rules.parse_count_matrix.output[0]
    output:
        cat = RES_DIR + "/tables/tcr_barcode.csv"
    script:
        "scripts/C_comb_barcodes_TCR.py"
        
rule augment_tcr_barcodes:
    """
    Merge HTO analysis with main data.
    Calculate binding concordance.
    Check clonotypes against VDJdb and IEDB data.
    """
    input:
        dat = rules.comb_barcodes_TCR.output.cat,
        vdj = WRK_DIR + "/tools/tcr_db.csv.gz",
        hto = RES_DIR + "/tables/hto.csv",
        gex = RES_DIR + "/tables/gex.txt"
    output:
        RES_DIR + "/tables/tcr_barcode.augmented.csv"
    shell:
        "python scripts/D_augment_tcr_barcodes.py \
            --data {input.dat} \
            --hto {input.hto} \
            --gex {input.gex} \
            --tcrdb {input.vdj} \
            --output {output}"

#################################################################
#            Perform grid search for optimal filters            #
#################################################################
rule eval_clonotypes:
    """
    Identify clonotypes with sufficient data to infer threshold based on.
    Plots barcode distribution for clonotypes with more than 10 GEMs.
    """
    input:
        data = rules.augment_tcr_barcodes.output[0],
        barcodes = EXPERIMENTAL_DESIGN_TEMPLATE
    params:
        plots = PLT_DIR + "/eval_clonotypes/%s/%d.pdf"
    output:
        data = RES_DIR + "/tables/tcr_barcode.valid.csv",
        done = touch(expand(PLT_DIR + "/eval_clonotypes/{flag}/dir.done", flag=['significant_match','significant_mismatch','insignificant']))
    conda:
        "envs/basic_dependencies.yaml"
    shell:
        "python scripts/F_comp_cred_specificities.py \
            --input {input.data} \
            --barcodes {input.barcodes} \
            --output {output.data} \
            --plots {params.plots}"


rule grid_search:
    """
    Run grid search on UMI values to define optimal set of thresholds.
    """
    input:
        valid_df = rules.eval_clonotypes.output.data
    output:
        grid = RES_DIR + "/eval_clonotypes/grid_search/{ext_thr}.csv"
    shell:
        "python scripts/G1_grid_search.py \
            --input {input.valid_df} \
            --ext_thr {wildcards.ext_thr} \
            --output {output.grid}"


rule extract_optimal_threshold:
    """
    Plot grid and extract optimal thresholds.
    """
    input:
        valid = rules.eval_clonotypes.output.data,
        grids = expand(rules.grid_search.output.grid, ext_thr=[0,1,2])
    output:
        plots = expand(PLT_DIR + "/eval_clonotypes/grid_search/grid.{ext}", ext=["pdf", "png"]),
        opt_thr = RES_DIR + "/eval_clonotypes/threshold/opt.csv"
    shell:
        "python scripts/G2_extract_optimal_threshold.py \
            --data {input.valid} \
            --grids {input.grids} \
            --output {output.opt_thr} \
            --plot {output.plots}"


#################################################################
#                             Plots                             #
#################################################################
rule get_filters:
    input:
        opt_thr = rules.extract_optimal_threshold.output.opt_thr,
        valid = rules.eval_clonotypes.output.data
    params:
        flt = '{filtering_set}'
    output:
        lbl = RES_DIR + "/eval_clonotypes/threshold/{filtering_set}.yaml",
        flt = RES_DIR + "/eval_clonotypes/threshold/{filtering_set}.csv"
    shell:
        "python scripts/H_set_filters.py \
            --data {input.valid} \
            --opt-thr {input.opt_thr} \
            --setting {params.flt} \
            --labels {output.lbl} \
            --filters {output.flt}"
        
rule filter_impact_staircase:
    input:
        df = rules.eval_clonotypes.output.data,
        lbl = rules.get_filters.output.lbl,
        flt = rules.get_filters.output.flt
    output:
        png = PLT_DIR + "/specificity_matrix/{filtering_set}/total.png"
    params:
        lambda wildcards, output: os.path.dirname(output.png)
    conda:
        "envs/basic_dependencies.yaml"
    shell:
        "python scripts/I_filter_impact.staircase.py \
            --data {input.df} \
            --labels {input.lbl} \
            --filters {input.flt} \
            --out-dir {params}"
        
rule ok:
    input:
        fn = expand(rules.filter_impact_staircase.output.png, filtering_set=['indv','comb'])
    output:
        tag = touch(EXP_DIR + '/done.ok')
    
