#################################################################
#                            Targets                            #
#################################################################
TARGET['plot_venn_MHC_AB'] = PLT_DIRECTORY + "/venn/MHC_AB.pdf"
TARGET['plot_venn_TCR_MHC_AB'] = expand(PLT_DIRECTORY + "/venn/TCR_MHC_AB.{sorting}.png", sorting=sorting_set)
#TARGET['plot_venn_HLA_concordance'] = PLT_DIRECTORY + "/venn/HLA_concordance.pdf"
TARGET['plot_confusion_multiplets'] = expand(PLT_DIRECTORY + "/confusion_multiplets/{sorting}.png", sorting=sorting_set)
#TARGET['specificity_matrices'] = expand(PLT_DIRECTORY + "/specificity_matrix/{plot_type}/{min_bc}.{min_tcr}.{ecs}.{ess}.pdf",
#                                        plot_type=config['plot_type'], min_bc=config['bc_thresh'],
#                                        min_tcr=config['tcr_thresh'],
#                                        ecs='ecs_False', ess=['ess_False'])

#################################################################
#                             Rules                             #
#################################################################
rule plot_venn_MHC_AB:
    input:
        MHC_DIRECTORY + "/count/brc.augmented.csv" #MHC_DIRECTORY + "/augmented/brc.csv"
    output:
        PLT_DIRECTORY + "/venn/MHC_AB.pdf"
    conda:
        "../envs/venn.yaml"
    script:
        "../../scripts/plot_venn_MHC_AB.py"

rule plot_venn_TCR_MHC_AB:
    input:
        MHC_DIRECTORY + "/count/brc.augmented.csv", #MHC_DIRECTORY + "/augmented/brc.csv",
        TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv"
    output:
        w = PLT_DIRECTORY + "/venn/TCR_MHC_AB.{sorting}.png",
        u = PLT_DIRECTORY + "/venn/TCR_MHC_AB.{sorting}.unweighted.png"
    params:
        label = "{sorting}"
    conda:
        "../envs/venn.yaml"
    priority:
        100
    script:
        "../../scripts/plot_venn_TCR_MHC_AB.py"

rule plot_venn_HLA_concordance:
	input:
		CAT_DIRECTORY + "/tables/tcr_barcode.csv",
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.imputed.csv"
	output:
		PLT_DIRECTORY + "/venn/HLA_concordance.pdf"
	conda:
		"../envs/venn.yaml"
	script:
		"../../scripts/plot_venn_HLA_concordance.py"

rule plot_confusion_multiplets:
    input:
        CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv" # Not imputed?
    output:
        PLT_DIRECTORY + "/confusion_multiplets/{sorting}.png" #{min_bc}
    params:
        sorting = '{sorting}'
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/plot_confusion_multiplets.py"

#rule plot_specificity_matrix:
#    input:
#        CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv" #.imputed
#    output:
#        plots = expand(PLT_DIRECTORY + "/specificity_matrix/{plot_type}/{{min_bc}}.{{min_tcr}}.{{ecs}}.{{ess}}.pdf", plot_type=config['plot_type'])
#    params:
#        plot_type = config['plot_type'],
#        min_brc = "{min_bc}",
#        min_tcr = "{min_tcr}",
#        ecs = "{ecs}",
#        ess = "{ess}"
#    conda:
#        "../envs/basic_dependencies.yaml"
#    script:
#        "../../scripts/E_plot_specificity_matrix.py"