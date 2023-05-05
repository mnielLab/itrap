
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

RAW_DATA = "experiments/exp13/run1/data/raw.csv.gz" # OBS! Change this to fit your path

WRK_DIR = workflow.basedir
EXP_DIR = os.path.join(WRK_DIR, "experiments", config["exp"], config["run"])
RES_DIR = os.path.join(EXP_DIR, "res")
PLT_DIR = os.path.join(EXP_DIR, "plt")


rule all:
    input: EXP_DIR + '/done.ok'


#################################################################
#            Perform grid search for optimal filters            #
#################################################################
rule eval_clonotypes:
    """
    Identify clonotypes with sufficient data to infer threshold based on.
    Plots barcode distribution for clonotypes with more than 10 GEMs.
    """
    input:
        data = RAW_DATA
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
    
