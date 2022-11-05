import pandas as pd
import os
#################################################################
#                            Targets                            #
#################################################################
if config['demultiplexed'] == False:
    TARGET['mkfastq'] = TCR_DIRECTORY + '/fastq/Stats/DemuxSummaryF1L1.txt' # OBS! Check if this path makes sense!

TARGET['ranger'] = expand(TCR_DIRECTORY + '/cellranger_{sorting}.done', sorting=sorting_set), # I expect this to be necessary
#TARGET['gex_filtering'] = CAT_DIRECTORY + "/tables/gex_filtering.txt"

#################################################################
#                             Rules                             #
#################################################################
rule mk_cellranger_configs:
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

rule clean_augment_tcr:
    input:
        expand(TCR_DIRECTORY + '/cellranger_{sorting}.done', sorting=sorting_set)
    params:
        contig = expand(TCR_DIRECTORY + '/cellranger_{sorting}/outs/multi/vdj_t/all_contig_annotations.csv', sorting=sorting_set),
        clonot = TCR_DIRECTORY + f'/cellranger_{total}/outs/per_sample_outs/cellranger_{total}/vdj_t/consensus_annotations.csv'
    output:
        TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
        TCR_DIRECTORY + "/reports/gems/gem_counts.json"
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/00_clean_augment_tcr.py"

#rule link_cellranger_contig_annotations_output:
#	input:
#		expand(TCR_DIRECTORY + '/cellranger_{sorting}.done', sorting=sorting.keys())
#	params:
#		*expand(TCR_DIRECTORY + '/cellranger_{sorting}/outs/multi/vdj_t/all_contig_annotations.csv', sorting=sorting.keys()) # Unpacking the list
#	output:                                                                                                                                                                      
#		TCR_DIRECTORY + '/all_contig_annotations.csv'
#	shell:
#		"""
#		./workflow/wrappers/cat_name_csv.sh {params} {output}
#		"""
#
#rule link_cellranger_clonotypes_output:
#	input:
#		expand(TCR_DIRECTORY + '/cellranger_{sorting}.done', sorting=sorting.keys())
#	params:
#		*expand(TCR_DIRECTORY + '/cellranger_{sorting}/outs/per_sample_outs/cellranger_{sorting}/vdj_t/consensus_annotations.csv', sorting=sorting.keys())                                           
#	output:                                                                                                                                                                      
#		TARGET['clonot']
#	run:
#		if len(input) > 1:
#			shell("cat {params} > {output}")
#		else:
#			shell("ln -sr {params} {output}")

#################################################################
#                        Gene expression                        #
#################################################################
ruleorder: filter_gex_data > no_gex_data

rule filter_gex_data:
    input:
        TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/count/raw_feature_bc_matrix.h5'
    output:
        CAT_DIRECTORY + "/eval_clonotypes/threshold/gex.txt",
        PLT_DIRECTORY + "/GEX/filtering/violin.png",
        PLT_DIRECTORY + "/GEX/filtering/scatter.png"
    conda:
        "../envs/seurat_gex_filtering.yaml"
    shell:
        "Rscript ./scripts/seurat_GEX_filtering.R {input} {output}"
        
rule no_gex_data:
    output:
        touch(CAT_DIRECTORY + "/eval_clonotypes/threshold/gex.txt")

#################################################################
#                Barcodes analyzed by Cellranger                #
#################################################################
rule link_cellranger_barcodes_output:
    input:
        TCR_DIRECTORY + f'/cellranger_{total}.done'
    params:
        brc = TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz',
        gem = TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz',
        mtx = TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz'
    output:
        brc = MHC_DIRECTORY + "/count/features.10x.tsv.gz",
        gem = MHC_DIRECTORY + "/count/barcodes.10x.tsv.gz",
        mtx = MHC_DIRECTORY + "/count/matrix.10x.mtx.gz",
        #tmp = touch(MHC_DIRECTORY + "/count/link.done")
    run:
        if os.path.isfile(params.brc) & os.path.isfile(params.gem) & os.path.isfile(params.mtx):
            shell("ln -sr {params.brc} {output.brc}")
            shell("ln -sr {params.gem} {output.gem}")
            shell("ln -sr {params.mtx} {output.mtx}")
            
rule link_cellranger_gems:
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





