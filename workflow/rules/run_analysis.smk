import numpy as np
import glob
#################################################################
#                            Targets                            #
#################################################################

TARGET['plot_thr_specificity_matrix'] = expand(PLT_DIRECTORY + "/specificity_matrix/peptide_per_clonotype_by_gem_size/{filtering_set}/plots.done",
                                               filtering_set=['indv','comb'])

TARGET['plot_thr_similarity'] = expand(PLT_DIRECTORY + '/eval_clonotypes/similarity/{filtering_set}/{plot_type}.png', #"/eval_clonotypes/similarity/{thrs_vals}/{plot_type}.pdf",
                                       filtering_set=['indv','comb'],
                                       plot_type=['score_per_pep','score_pooled','auc','roc'])

TARGET['plot_thr_impact'] = expand(PLT_DIRECTORY + "/eval_filters/{filtering_set}/sco.png", filtering_set=['indv','comb'])

#################################################################
#                        GLOBAL variables                       #
#################################################################
# OBS! This range is heuristic. May require changes.
thr_conv_dct = {str(k):v for k,v in enumerate(2**np.linspace(-0.4,3,50))} # delta umi mhc

# Number of times to sample TCRs for similarity scores.
random_sample=1 #np.arange(30)

#################################################################
#                             Rules                             #
#################################################################

rule eval_clonotypes:
    """
    Identify clonotypes with sufficient data to infer threshold based on.
    Plots barcode distribution for clonotypes with more than 10 GEMs.
    """
    input:
        data = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
        barcodes = LIB_DIRECTORY + '/barcode_specificity_annotations.xlsx'
    params:
        plots = PLT_DIRECTORY + "/eval_clonotypes/%s/%d.pdf"
    output:
        data = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        done = touch(expand(PLT_DIRECTORY + "/eval_clonotypes/{flag}/dir.done",flag=['significant_match','significant_mismatch','insignificant']))
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/F_comp_cred_specificities.py"


rule grid_search:
    """
    Run grid search on UMI values to define optimal set of thresholds.
    """
    input:
        valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
    params:
        ext_thr = lambda wildcards: thr_conv_dct[wildcards.external_threshold]
    output:
        grid = CAT_DIRECTORY + "/eval_clonotypes/grid_search/{external_threshold}.csv"
    script:
        "../../scripts/G1_grid_search.py"


rule extract_optimal_threshold:
    """
    Plot grid and extract optimal thresholds.
    """
    input:
        valid = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        grids = expand(CAT_DIRECTORY + "/eval_clonotypes/grid_search/{external_threshold}.csv", external_threshold=thr_conv_dct.keys())
    output:
        plots = expand(PLT_DIRECTORY + "/eval_clonotypes/grid_search/opt_thr/grid.{ext}", ext=["pdf", "png"]),
        opt_thr = CAT_DIRECTORY + "/eval_clonotypes/threshold/opt.csv"
    script:
        "../../scripts/G2_extract_optimal_threshold.py"


rule plot_thr_specificity_matrix:
    """
    Plot specificity: pMHC per clonotype.
    """
    input:
        thresholds = CAT_DIRECTORY + "/eval_clonotypes/threshold/{thrs_vals}.csv",
        original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
        valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
    output:
        plots = PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{thrs_vals}/{thrs_levl}/{plot_type}.{ext}"
    params:
        plot_type = lambda wildcards: wildcards.plot_type,
        thrs_levl = lambda wildcards: wildcards.thrs_levl
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/E_plot_specificity_matrix-Copy1.thr.py"
        
rule plot_sample_specificity_matrix:
    """
    Plot specificity per sample given by cell hashing.
    """
    input:
        thresholds = CAT_DIRECTORY + "/eval_clonotypes/threshold/{thrs_vals}.csv",
        original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
        valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
    output:
        #plots = PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{thrs_vals}/{plot_type}/{data_type}.{thrs_type}.{ext}"
        touch(PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{thrs_vals}/{plot_type}_{data_type}_{thrs_type}/samples.done")
    params:
        plot_type = lambda wildcards: wildcards.plot_type,
        data_type = lambda wildcards: wildcards.data_type,
        thrs_type = lambda wildcards: wildcards.thrs_type
    conda:
        "../envs/basic_dependencies.yaml"
    #priority:
    #    90
    script:
        "../../scripts/E_plot_specificity_matrix.sample.py"
    

#rule kernel_similarity_between_all_cdr3s:
#    input:
#        s2s = config['s2s'], #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel",
#        blf = config['blf'], #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50",
#        qij = config['qij'], #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij",
#        tmp = TCR_DIRECTORY + f'/cellranger_{total}.done'
#    params:
#        df = TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/vdj_t/all_contig_annotations.csv'
#    output:
#        A = CAT_DIRECTORY + "/similarity_assessment/cdr3_a.csv",
#        B = CAT_DIRECTORY + "/similarity_assessment/cdr3_b.csv"
#    script:
#        "../../scripts/kernel_similarity_between_all_cdr3s.py"

rule get_all_cdr3_seqs:
    """
    Extract all unique CDR3 sequences for similarity computations.
    """
    input:
        TCR_DIRECTORY + f'/cellranger_{total}.done'
    params:
        TCR_DIRECTORY + f'/cellranger_{total}/outs/multi/vdj_t/all_contig_annotations.csv',
        chain = '{chain}'
    output:
        CAT_DIRECTORY + "/similarity_assessment/{chain}_cdr3_all.txt"
    script:
        "../../scripts/K_get_all_cdr3_seqs.py"

# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint split_cdr3_seqs_in_chuncks:
    """
    For parallelization, split list of unique CDR3 sequences into smaller files.
    """
    input:
        CAT_DIRECTORY + "/similarity_assessment/{chain}_cdr3_all.txt"
    params:
        chain = CAT_DIRECTORY + "/similarity_assessment/{chain}"
    output:
        directory(CAT_DIRECTORY + "/similarity_assessment/chunks/{chain}/")
    run:
        lst = np.loadtxt(input[0], dtype='str')
        # Write chunck files for parallel processing
        n = 200
        for i in range(0, len(lst), n):
            np.savetxt(output[0] + str(i) + '.txt', lst[i:i + n], fmt='%s')
            
rule kernel_similarity_between_all_cdr3s:
    """
    Compute similarity using TCRmatch kernel method.
    """
    priority: 100
    input:
        s2s = config['s2s'],
        blf = config['blf'],
        qij = config['qij'],
        chu = CAT_DIRECTORY + "/similarity_assessment/chunks/{chain}/{chunk}.txt",
        cdr = CAT_DIRECTORY + "/similarity_assessment/{chain}_cdr3_all.txt"
    params:
        chunk_i = '{chunk}',
        chain = '{chain}'
    output:
        chunk = CAT_DIRECTORY + "/similarity_assessment/chunks/{chain}/{chunk}.csv"
    script:
        "../../scripts/K_kernel_similarity_between_all_cdr3s.py"
        
def aggregate_similarity_computations(wildcards):
    checkpoint_output = checkpoints.split_cdr3_seqs_in_chuncks.get(**wildcards).output[0]
    chunks = glob_wildcards(os.path.join(checkpoint_output, "{chunk}.txt")).chunk
    return expand(CAT_DIRECTORY + "/similarity_assessment/chunks/{chain}/{chunk}.csv", chain=wildcards.chain, chunk=chunks)
        
rule aggregate_similarity_computations:
    input:
        aggregate_similarity_computations
    output:
        CAT_DIRECTORY + "/similarity_assessment/cdr3_{chain}.csv"
    shell:
        "cat {input} > {output}"
        

checkpoint index_similarity_per_peptide:
    """
    Extract index of clonotypes for each peptide to do similarity computations.
    """
    input:
        valid = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        idx_tra = CAT_DIRECTORY + "/similarity_assessment/a_cdr3_all.txt",
        idx_trb = CAT_DIRECTORY + "/similarity_assessment/b_cdr3_all.txt"
    output:
        directory(CAT_DIRECTORY + "/similarity_assessment/tmp_cmd/")
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/G_idx_sim_per_pep.py"
        
        
rule split_similarity_per_peptide:
    """
    Extract similarity scores by index for each pMHC and each chain (alpha & beta).
    """
    input:
        idx = CAT_DIRECTORY + "/similarity_assessment/tmp_cmd/{chain}.{pep}.txt",
        sim = CAT_DIRECTORY + "/similarity_assessment/cdr3_{chain}.csv"
    params:
        tmp = lambda wildcards, input: input[0].strip('txt') + 'idx'
    output:
        CAT_DIRECTORY + "/similarity_assessment/tmp_csv/{chain}.{pep}.csv"
    shell:
        """
        while read rownum; do
            printf '%.12d\n' "$rownum"
        done <{input.idx} >{params.tmp}
        
        join -1 2 -2 1 <(nl {params.tmp} | sort -k 2,2) <(nl -w 12 -n rz {input.sim}) | sort -k 2,2n | cut -d ' ' -f 3- > {output}
        """
        
rule compute_similarity_per_peptide:
    """
    Compute intra and inter similarity per peptide for different sets of filters.
    Computes for filters combined or individually (filtering set).
    Allows for multiple replicas in case of a small dataset (random sample).
    """
    input:
        opt_thr = CAT_DIRECTORY + "/eval_clonotypes/threshold/opt.csv",
        valid = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        gex = CAT_DIRECTORY + "/eval_clonotypes/threshold/gex.txt",
        sim_tra = CAT_DIRECTORY + "/similarity_assessment/tmp_csv/a.{pep}.csv",
        sim_trb = CAT_DIRECTORY + "/similarity_assessment/tmp_csv/b.{pep}.csv",
        idx_tra_sub = CAT_DIRECTORY + "/similarity_assessment/tmp_cmd/a.{pep}.txt",
        idx_trb_sub = CAT_DIRECTORY + "/similarity_assessment/tmp_cmd/b.{pep}.txt",
        idx_tra = CAT_DIRECTORY + "/similarity_assessment/a_cdr3_all.txt",
        idx_trb = CAT_DIRECTORY + "/similarity_assessment/b_cdr3_all.txt"
    params:
        flt = '{filtering_set}',
        rnd = '{random_sample}',
        pep = '{pep}'
    output:
        sim = CAT_DIRECTORY + "/similarity_assessment/sim_pep/{filtering_set}/{random_sample}/{pep}.csv"
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/G_similarity_per_peptide.py"
        
def aggregate_similarity_per_pep(wildcards):
    checkpoint_output = checkpoints.index_similarity_per_peptide.get(**wildcards).output[0]
    pep = glob_wildcards(os.path.join(checkpoint_output, "{chain}.{pep}.txt")).pep
    return expand(CAT_DIRECTORY + "/similarity_assessment/sim_pep/{filtering_set}/{random_sample}/{pep}.csv",
                  filtering_set=wildcards.filtering_set, random_sample=wildcards.random_sample, pep=pep)
        
rule agg_similarity_per_peptide:
    """
    Combine all similarity scores across all peptides still partitioned by sampling if chosen.
    """
    input:
        aggregate_similarity_per_pep
    params:
        flt = '{filtering_set}',
        rnd = '{random_sample}'
    output:
        CAT_DIRECTORY + "/similarity_assessment/plt_df/{filtering_set}/{random_sample}.csv"
    shell:
        "cat {input} > {output}"
        
#rule test_thr_similarity:
#    input:
#        opt_thr = CAT_DIRECTORY + "/eval_clonotypes/threshold/opt.csv", #opt {thrs_vals}
#        valid = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
#        gex = CAT_DIRECTORY + "/eval_clonotypes/threshold/gex.txt",
#        sim_tra = CAT_DIRECTORY + "/similarity_assessment/cdr3_a.csv",
#        sim_trb = CAT_DIRECTORY + "/similarity_assessment/cdr3_b.csv",
#        idx_tra = CAT_DIRECTORY + "/similarity_assessment/a_cdr3_all.txt",
#        idx_trb = CAT_DIRECTORY + "/similarity_assessment/b_cdr3_all.txt"
#    params:
#        flt = '{filtering_set}',
#        rnd = '{random_sample}'
#    output:
#        sim = CAT_DIRECTORY + "/similarity_assessment/plt_df/{filtering_set}/{random_sample}.csv"
#        #auc = CAT_DIRECTORY + "/similarity_assessment/plt_df/auc.{filtering_set}.csv"
#    #priority:
#    #90
#    conda:
#        "../envs/basic_dependencies.yaml"
#    script:
#        "../../scripts/G_test_inter_v_inter_similarity.py"
        
rule comb_thr_similarity:
    input:
        expand(CAT_DIRECTORY + "/similarity_assessment/plt_df/{{filtering_set}}/{random_sample}.csv", random_sample=random_sample)
    output:
        CAT_DIRECTORY + "/similarity_assessment/plt_df/{filtering_set}.sim.csv"
    shell:
        """
        cat {input} > {output}
        """
        
rule comp_auc_similarity:
    input:
        sim = CAT_DIRECTORY + "/similarity_assessment/plt_df/{filtering_set}.sim.csv"
    output:
        auc = CAT_DIRECTORY + "/similarity_assessment/plt_df/{filtering_set}.auc.csv"
    script:
        "../../scripts/G_comp_inter_v_inter_auc.py"
        
rule plot_thr_similarity:
    input:
        sim = CAT_DIRECTORY + "/similarity_assessment/plt_df/{filtering_set}.sim.csv",
        auc = CAT_DIRECTORY + "/similarity_assessment/plt_df/{filtering_set}.auc.csv"
    params:
        flt = '{filtering_set}'
    output:
        score_per_pep = PLT_DIRECTORY + "/eval_clonotypes/similarity/{filtering_set}/score_per_pep.png",
        score_pooled = PLT_DIRECTORY + "/eval_clonotypes/similarity/{filtering_set}/score_pooled.png",
        auc = PLT_DIRECTORY + "/eval_clonotypes/similarity/{filtering_set}/auc.png",
        roc = PLT_DIRECTORY + "/eval_clonotypes/similarity/{filtering_set}/roc.png",
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/G_plot_inter_v_inter_similarity.py"
        
rule get_filters:
    input:
        opt_thr = CAT_DIRECTORY + "/eval_clonotypes/threshold/opt.csv",
        valid = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
    params:
        flt = '{filtering_set}'
    output:
        lbl = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.yaml",
        flt = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.csv"
    script:
        "../../scripts/H_set_filters.py"
        
rule filter_impact_staircase:
    input:
        df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        lbl = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.yaml",
        flt = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.csv"
    output:
        touch(PLT_DIRECTORY + "/specificity_matrix/peptide_per_clonotype_by_gem_size/{filtering_set}/plots.done")
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/I_filter_impact.staircase.py"
        
rule filter_impact_bar:
    input:
        df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        flt = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.yaml",
        idx = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.csv",
        auc = CAT_DIRECTORY + "/similarity_assessment/plt_df/{filtering_set}.auc.csv"
    output:
        sco = PLT_DIRECTORY + "/eval_filters/{filtering_set}/sco.png",
        old = PLT_DIRECTORY + "/eval_filters/{filtering_set}/old.png",
        tcr = PLT_DIRECTORY + "/eval_filters/{filtering_set}/tcr.png",
        pep = PLT_DIRECTORY + "/eval_filters/{filtering_set}/pep.png"
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/filtering_impact.bar.py"
