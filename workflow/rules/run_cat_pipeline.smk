#################################################################
#                            Targets                            #
#################################################################
TARGET['tcr_brc_cln'] = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"

#################################################################
#                             Rules                             #
#################################################################

#################################################################
#           Merge Barcodes analyzed by both platforms           #
#################################################################
# How to compile output from KMA and cellranger?
# Use rule inheritance:
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-inheritance
# Use Handling Ambiguous Rules
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#handling-ambiguous-rules

ruleorder: parse_count_matrix > parse_count_matrix_2 > parse_count_matrix_3
    
rule parse_count_matrix:
    input:
        brc = expand(MHC_DIRECTORY + "/count/features.{ext}", ext=['kma.tsv','10x.tsv.gz']),
        gem = expand(MHC_DIRECTORY + "/count/barcodes.{ext}", ext=['kma.tsv','10x.tsv.gz']),
        mtx = expand(MHC_DIRECTORY + "/count/matrix.{ext}", ext=['kma.mtx','10x.mtx.gz']),
        cnt = [MHC_DIRECTORY + "/mapping/umi.tsv", MHC_DIRECTORY + "/count/gem_labels.10x.csv"],
        lbl = LIB_DIRECTORY + '/barcode_specificity_annotations.xlsx',
        ann = WRK_DIR + '/tools/detected_responses_annotation.xlsx'
    params:
        mhc_custom = mhc_custom,
        hsh_custom = hsh_custom,
        mrk_custom = mrk_custom
    output:
        MHC_DIRECTORY + "/count/brc.augmented.csv"
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/parse_count_matrix.py"
        
rule parse_count_matrix_2:
    input:
        brc = expand(MHC_DIRECTORY + "/count/features.{ext}", ext=['kma.tsv']),
        gem = expand(MHC_DIRECTORY + "/count/barcodes.{ext}", ext=['kma.tsv']),
        mtx = expand(MHC_DIRECTORY + "/count/matrix.{ext}", ext=['kma.mtx']),
        cnt = [MHC_DIRECTORY + "/mapping/umi.tsv"],
        lbl = LIB_DIRECTORY + '/barcode_specificity_annotations.xlsx',
        ann = WRK_DIR + '/tools/detected_responses_annotation.xlsx'
    params:
        mhc_custom = mhc_custom,
        hsh_custom = hsh_custom,
        mrk_custom = mrk_custom
    output:
        MHC_DIRECTORY + "/count/brc.augmented.csv"
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/parse_count_matrix.py"
        
rule parse_count_matrix_3:
    input:
        brc = expand(MHC_DIRECTORY + "/count/features.{ext}", ext=['10x.tsv.gz']),
        gem = expand(MHC_DIRECTORY + "/count/barcodes.{ext}", ext=['10x.tsv.gz']),
        mtx = expand(MHC_DIRECTORY + "/count/matrix.{ext}", ext=['10x.mtx.gz']),
        cnt = [MHC_DIRECTORY + "/count/gem_labels.10x.csv"],
        lbl = LIB_DIRECTORY + '/barcode_specificity_annotations.xlsx',
        ann = WRK_DIR + '/tools/detected_responses_annotation.xlsx'
    params:
        mhc_custom = mhc_custom,
        hsh_custom = hsh_custom,
        mrk_custom = mrk_custom
    output:
        MHC_DIRECTORY + "/count/brc.augmented.csv"
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/parse_count_matrix.py"

#################################################################
#                     Merge Barcodes and TCR                    #
#################################################################

rule comb_barcodes_TCR: # Overflødig / eller erstat med run.
    input:
        TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
        MHC_DIRECTORY + "/count/brc.augmented.csv" #MHC_DIRECTORY + "/augmented/brc.csv"
    output:
        CAT_DIRECTORY + "/tables/tcr_barcode.csv"
    script:
        "../../scripts/C_comb_barcodes_TCR.py"
        
#rule prep_seurat_HTO_analysis:
#    input:
#        df = CAT_DIRECTORY + "/tables/tcr_barcode.csv",
#        hto = LIB_DIRECTORY + '/barcode_specificity_annotations.xlsx',
#        brc = expand(MHC_DIRECTORY + "/count/features.{ext}", ext=['10x.tsv.gz']),
#        gem = expand(MHC_DIRECTORY + "/count/barcodes.{ext}", ext=['10x.tsv.gz']),
#        mtx = expand(MHC_DIRECTORY + "/count/matrix.{ext}", ext=['10x.mtx.gz'])
#    output:
#        temp(CAT_DIRECTORY + "/tables/hto_count_matrix.rds")
#    conda:
#        "../envs/prep_seurat.yaml"
#    script:
#        "../../scripts/C1_prep_seurat.py"
        
rule prep_seurat_HTO_analysis:
    input:
        df = CAT_DIRECTORY + "/tables/tcr_barcode.csv",
        hto = LIB_DIRECTORY + '/barcode_specificity_annotations.xlsx',
        brc = MHC_DIRECTORY + "/count/features.10x.tsv.gz",
        gem = MHC_DIRECTORY + "/count/barcodes.10x.tsv.gz",
        mtx = MHC_DIRECTORY + "/count/matrix.10x.mtx.gz"
    params:
        CAT_DIRECTORY + "/tables/hto_prep.done"
    output:
        temp(CAT_DIRECTORY + "/tables/hto_count_matrix.rds")
    conda:
        "../envs/prep_seurat.yaml"
    script:
        "../../scripts/C1_prep_seurat-Copy1.py"

ruleorder: seurat_HTO_analysis > no_seurat_HTO_analysis
rule seurat_HTO_analysis:
    input:
        CAT_DIRECTORY + "/tables/hto_count_matrix.rds",
        CAT_DIRECTORY + "/tables/hto_prep.done"
    output:
        out_file = CAT_DIRECTORY + "/tables/hto.csv",
        ridge_plot = PLT_DIRECTORY + "/HTO/ridgeplot.png",
        violin_plot = PLT_DIRECTORY + "/HTO/violinplot.png",
        heatmap = PLT_DIRECTORY + "/HTO/heatmap.png"
    conda:
        "../envs/seurat_hto.yaml"
    shell:
        "Rscript ./scripts/seurat_HTO_analysis.R {input} {output.out_file} {output.ridge_plot} {output.violin_plot} {output.heatmap}"
 
rule no_seurat_HTO_analysis:
    output:
        out_file = CAT_DIRECTORY + "/tables/hto.csv"
    run:
        import pandas as pd
        
        df = pd.DataFrame(columns=['gem','seurat','umi_count_hto','feature_hto','hto_max_id','hto_sec_id',
                                   'hto_margin','hto_classification','hto_global_class','hash_id'])
        df.to_csv(output.out_file, index=False)

rule clean_tcr_barcodes:
    input:
        dat = CAT_DIRECTORY + "/tables/tcr_barcode.csv", #rules.comb_barcodes_TCR.output, #
        vdj = WRK_DIR + "/tools/VDJdb.csv",
        hto = CAT_DIRECTORY + "/tables/hto.csv"
    params:
        TCR_DIRECTORY + "/library/clone_rearrangement.tsv" # Change position!
    output:
        CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
        CAT_DIRECTORY + "/reports/gems/gem_counts.json"
    script:
        "../../scripts/D_clean_tcr_barcodes.py"
        







