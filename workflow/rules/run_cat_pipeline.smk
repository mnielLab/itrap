#################################################################
#                            Targets                            #
#################################################################
"""
Add intermediate targets to pipeline to test partial run, e.g. add:

TARGET['tcr_brc_cln'] = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"
"""

TARGET['tcr_brc_cln'] = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"

#################################################################
#                             Rules                             #
#################################################################

#################################################################
#           Merge Barcodes analyzed by both platforms           #
#################################################################

# Test the following methods to parse count matrix
ruleorder: parse_count_matrix > parse_count_matrix_2 > parse_count_matrix_3
    
rule parse_count_matrix:
    """
    Merge count matrices of custom barcode reads and barcode reads identified with 10x cellranger.
    Annotate types of barcode: pMHC, cell hashing, or surface markers.
    Annotate most abundant barcode per GEM and other summary metrics.
    """
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
    """
    Run this rule if barcode reads origin from custom barcodes and are mapped using KMA.
    Annotate types of barcode: pMHC, cell hashing, or surface markers.
    Annotate most abundant barcode per GEM and other summary metrics.
    """
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
    """
    Run this rule if barcode reads origin from 10x barcodes and are mapped using Cellranger.
    Annotate types of barcode: pMHC, cell hashing, or surface markers.
    Annotate most abundant barcode per GEM and other summary metrics.
    """
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

rule comb_barcodes_TCR:
    """
    Merge barcode and TCR data based on GEM barcode.
    """
    input:
        TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
        MHC_DIRECTORY + "/count/brc.augmented.csv"
    output:
        CAT_DIRECTORY + "/tables/tcr_barcode.csv"
    script:
        "../../scripts/C_comb_barcodes_TCR.py"
        
rule prep_seurat_HTO_analysis:
    """
    Convert count matrix to RDS format for Seurat analysis.
    """
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
    """
    Analyze cell hashing barcodes with Seurat HTO method.
    """
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
    """
    If no cell hashing barcodes are available, mimic the HTO output.
    """
    output:
        out_file = CAT_DIRECTORY + "/tables/hto.csv"
    run:
        import pandas as pd
        
        df = pd.DataFrame(columns=['gem','seurat','umi_count_hto','feature_hto','hto_max_id','hto_sec_id',
                                   'hto_margin','hto_classification','hto_global_class','hash_id'])
        df.to_csv(output.out_file, index=False)

rule clean_tcr_barcodes:
    """
    Merge HTO analysis with main data.
    Calculate binding concordance.
    Check clonotypes against VDJdb and IEDB data.
    """
    input:
        dat = CAT_DIRECTORY + "/tables/tcr_barcode.csv",
        vdj = WRK_DIR + "/tools/VDJdb.csv",
        hto = CAT_DIRECTORY + "/tables/hto.csv"
    params:
        TCR_DIRECTORY + "/library/clone_rearrangement.tsv"
    output:
        CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
        CAT_DIRECTORY + "/reports/gems/gem_counts.json"
    script:
        "../../scripts/D_clean_tcr_barcodes.py"
        
