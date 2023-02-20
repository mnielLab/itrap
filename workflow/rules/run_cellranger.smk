import pandas as pd
import os

#################################################################
#                            Targets                            #
#################################################################

"""
Add intermediate targets to pipeline to test partial run, e.g. add:

TARGET['tcr_clean'] = TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv"
"""

TARGET['tcr_clean'] = TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv"

#################################################################
#                             Rules                             #
#################################################################
rule mk_cellranger_configs:
    """
    Writes standard config files for running Cellranger multi:
    https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi#config
    """
    input:
        barcode_annotations = LIB_DIRECTORY + '/barcode_specificity_annotations.xlsx',
        barcode_library = config['commercial_brc'],
        brc_info = config['brc_info'], # load all_oligos.csv
        fq_files = LIB_DIRECTORY + "/fq_files.py"
    output:
        multi_config = TCR_DIRECTORY + '/config/{sorting}_multi.csv'
    params:
        feature_reference = TCR_DIRECTORY + '/config/feature.csv'
    script:
        "../../scripts/write_cellranger_multi_config_files.py"

rule run_cellranger:
    """
    Running Cellranger multi.
    Cellranger generates a system of files and directories from CWD.
    As both Snakemake and Cellranger generates the output directories, the pipeline breaks.
    The work-around is to generate a faux output (.done).
    """
    input:
        multi_config = TCR_DIRECTORY + '/config/{sorting}_multi.csv'
    params:
        rundir = TCR_DIRECTORY,
        cellranger = config['cellranger']
    output:
        touch(TCR_DIRECTORY + '/cellranger_{sorting}.done')
    message: "Running Cellranger multi {wildcards.sorting}"
    shell:
        """
        cd {params.rundir}
        
        {params.cellranger} multi --id cellranger_{wildcards.sorting} --csv {input.multi_config}
        """
        
rule label_tcr:
    input:
        expand(TCR_DIRECTORY + '/cellranger_{sorting}.done', sorting=sorting_set)
    params:
        contig = expand(TCR_DIRECTORY + '/cellranger_{sorting}/outs/multi/vdj_t/all_contig_annotations.csv', sorting=sorting_set),
        clonot = TCR_DIRECTORY + f'/cellranger_{total}/outs/per_sample_outs/cellranger_{total}/vdj_t/consensus_annotations.csv'
    output:
        TCR_DIRECTORY + "/augmented/tcr.labeled.csv"
    run:
        for filename in params.contig:
            if 'tot' in filename:
                INPUT = filename
            elif 'pos' in filename:
                POSITIVE_GEMS = filename
            elif 'neg' in filename:
                NEGATIVE_GEMS = filename
                
        # Load data
        df = pd.read_csv(INPUT)
        try:
            pos_df = pd.read_csv(POSITIVE_GEMS)
            pos_gems = set(pos_df.barcode)
        except:
            pos_gems = set(df.barcode)
        try:
            neg_df = pd.read_csv(NEGATIVE_GEMS)
            neg_gems = set(neg_df.barcode) - pos_gems
        except:
            neg_gems = set()
            
        # Label data
        df.loc[df.barcode.isin(pos_gems), 'label'] = 1
        df.loc[df.barcode.isin(neg_gems), 'label'] = 0
        
        df.to_csv(output[0], index=False)

rule clean_augment_tcr:
    """
    Reannotating clonotypes:
    - Clonotypes of identical TCRab AA sequences are merged
    - GEMs with no clonotypes may be assigned an clonotype in two ways:
      1. if the TCRab pairs uniquely match an existing clonotype 
      2. if the TCRab are represent a novel clonotype
    """
    input:
        contig = TCR_DIRECTORY + "/augmented/tcr.labeled.csv"
    params:
        clonot = TCR_DIRECTORY + f'/cellranger_{total}/outs/per_sample_outs/cellranger_{total}/vdj_t/consensus_annotations.csv'
    output:
        output = TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
        report = TCR_DIRECTORY + "/reports/gems/gem_counts.json"
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/A_clean_augment_tcr.py"


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
        TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/count/raw_feature_bc_matrix.h5'
    output:
        CAT_DIRECTORY + "/eval_clonotypes/threshold/gex.txt",
        PLT_DIRECTORY + "/GEX/filtering/violin.png",
        PLT_DIRECTORY + "/GEX/filtering/scatter.png"
    conda:
        "../envs/seurat_gex_filtering.yaml"
    shell:
        "Rscript ./scripts/B0_seurat_gex_analysis.R {input} {output}"
        
rule no_gex_data:
    """
    Generate a dummy GEX filtering file.
    """
    output:
        touch(CAT_DIRECTORY + "/eval_clonotypes/threshold/gex.txt")

#################################################################
#                Barcodes analyzed by Cellranger                #
#################################################################
rule link_cellranger_barcodes_output:
    """
    Symlink the multimer DNA-barcode count matrix.
    """
    input:
        TCR_DIRECTORY + f'/cellranger_{total}.done'
    params:
        brc = TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz',
        gem = TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz',
        mtx = TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz'
    output:
        brc = MHC_DIRECTORY + "/count/features.10x.tsv.gz",
        gem = MHC_DIRECTORY + "/count/barcodes.10x.tsv.gz",
        mtx = MHC_DIRECTORY + "/count/matrix.10x.mtx.gz"
    run:
        if os.path.isfile(params.brc) & os.path.isfile(params.gem) & os.path.isfile(params.mtx):
            shell("ln -sr {params.brc} {output.brc}")
            shell("ln -sr {params.gem} {output.gem}")
            shell("ln -sr {params.mtx} {output.mtx}")
            
rule link_cellranger_gems:
    """
    Label GEMs based on multimer sorting.
    - Multimer positive cells get label 1
    - Multimer negative cells get label 0
    - If no negatively sorted cells are included, all GEMs are labelled 1.
    """
    input:
        expand(TCR_DIRECTORY + '/cellranger_{sorting}.done', sorting=sorting.keys())
    params:
        filenames = expand(TCR_DIRECTORY + '/cellranger_{sorting}/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz', sorting=sorting.keys())
    output:
        gems = MHC_DIRECTORY + "/count/gem_labels.10x.csv"
    run:
        lst = list()
        for filename in params.filenames:
            if os.path.isfile(filename):
                if 'pos' in filename:
                    label = 1
                if 'tot' in filename:
                    label = 1
                if 'neg' in filename:
                    label = 0

                df = pd.read_csv(filename, sep='\t', header=None, names=['gem'])
                df['tmp1'] = np.nan
                df['tmp2'] = np.nan
                df['label'] = label
                lst.append(df)
            
        tmp = pd.concat(lst)
        tmp.to_csv(output.gems, index=False, header=False)
        
        #assert tmp.duplicated().any() == False