rule clean_augment_tcr:
    input:
        TCR_DIRECTORY + '/cellranger_{sorting}.done'
    params:
        contig = TCR_DIRECTORY + '/cellranger_{sorting}/outs/multi/vdj_t/all_contig_annotations.csv', sorting=sorting.keys(),
        clonot = TCR_DIRECTORY + '/cellranger_{sorting}/outs/per_sample_outs/cellranger_{sorting}/vdj_t/consensus_annotations.csv'
    output:
        contig = TCR_DIRECTORY + "/augmented/tcr.clean.augmented.{sorting}.csv",
        report = TCR_DIRECTORY + "/reports/gems/gem_counts.json"
    script:
        "../../scripts/00_clean_augment_tcr.py"

rule link_cellranger_contig_annotations_output:
    input:
        *expand(TCR_DIRECTORY + "/augmented/tcr.clean.augmented.{sorting}.csv", sorting=sorting.keys()) # Unpacking the list
    output:
        TCR_DIRECTORY + '/all_contig_annotations.csv'
    shell:
        """
        ./workflow/wrappers/concat_cellranger_output.sh {input} {output}
        """

rule link_cellranger_clonotypes_output:
    input:
        expand(TCR_DIRECTORY + '/cellranger_{sorting}.done', sorting=sorting.keys())
    params:
        *expand(TCR_DIRECTORY + '/cellranger_{sorting}/outs/per_sample_outs/cellranger_{sorting}/vdj_t/consensus_annotations.csv', sorting=sorting.keys())                                           
    output:
        TARGET['clonot']
    run:
        if len(input) > 1:
            shell("cat {params} > {output}")
        else:
            shell("ln -sr {params} {output}")


