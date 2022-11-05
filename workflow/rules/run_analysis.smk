import numpy as np
import glob
#################################################################
#                            Targets                            #
#################################################################
#TARGET['grid_search'] = expand(CAT_DIRECTORY + "/eval_clonotypes/grid_search.{umi_count_mhc_rel}.csv",
#                               umi_count_mhc_rel=[0]) # list((np.linspace(0,0.2,20) * 10000).astype(int).astype(str))

# Testing
#TARGET['plot_thr_specificity_matrix'] = expand(PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{thrs_vals}/{plot_type}.{data_type}.{thrs_type}.{ext}",
#                                               thrs_vals=config['thrs_vals'],
#                                               plot_type=config['plot_type'],
#                                               data_type=['total'], #'train','test',
#                                               thrs_type=['raw','thr'],
#                                               ext=['pdf','png'])

#TARGET['plot_thr_specificity_matrix'] = expand(PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{thrs_vals}/{thrs_levl}/{plot_type}.{ext}",
#                                               thrs_vals=['opt'], #config['thrs_vals'],
#                                               thrs_levl=['1_thr_only','2_matching_hashing_only','3_specificity_multiplets_only','4_complete_tcrs_only'],
#                                               plot_type=config['plot_type'],
#                                               ext=['pdf','png'])
#


#TARGET['plot_sample_specificity'] = expand(PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{thrs_vals}/{plot_type}_{data_type}_{thrs_type}/samples.done",
#                                          thrs_vals=config['thrs_vals'],
#                                          plot_type=config['plot_type'],
#                                          data_type=['total'], #'train','test',
#                                          thrs_type=['raw','thr'])

TARGET['plot_thr_specificity_matrix'] = expand(PLT_DIRECTORY + "/specificity_matrix/peptide_per_clonotype_by_gem_size/{filtering_set}/plots.done",
                                               filtering_set=['indv','comb'])

TARGET['plot_thr_similarity'] = expand(PLT_DIRECTORY + '/eval_clonotypes/similarity/{filtering_set}/{plot_type}.png', #"/eval_clonotypes/similarity/{thrs_vals}/{plot_type}.pdf",
                                       filtering_set=['indv','comb'],
                                       plot_type=['score_per_pep','score_pooled','auc','roc'])

TARGET['plot_thr_impact'] = expand(PLT_DIRECTORY + "/eval_filters/{filtering_set}/sco.png", filtering_set=['indv','comb'])

#################################################################
#                             Rules                             #
#################################################################

rule eval_clonotypes: # Doesnt plot?
    input:
        CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
        LIB_DIRECTORY + '/barcode_specificity_annotations.xlsx'
    params:
        PLT_DIRECTORY + "/eval_clonotypes/%s/%d.pdf"
    output:
        CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        touch(expand(PLT_DIRECTORY + "/eval_clonotypes/{flag}/dir.done",flag=['significant_match','significant_mismatch','insignificant']))
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/F_comp_cred_specificities.py"
        
rule get_minimum_threshold:
    input:
        data = CAT_DIRECTORY + "/tables/tcr_barcode.csv"
    output:
        thrs = CAT_DIRECTORY + "/eval_clonotypes/threshold_min.csv",
        plot = expand(PLT_DIRECTORY + "/eval_clonotypes/opt_thr/min_thr.{tp}.pdf", tp=['brc','tcr'])
    script:
        "../../scripts/plot_bc_umi_dist_against_bg.py"
        
rule get_parameter_distributions:
    input:
        data = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"
    params:
        plot = expand(PLT_DIRECTORY + "/eval_clonotypes/param_dist/%s.{ext}", ext=["pdf", "png"])
    output:
        touch(PLT_DIRECTORY + "/eval_clonotypes/param_dist/param_dist.done"),
        prms = CAT_DIRECTORY + "/eval_clonotypes/param_dist.csv"
    conda:
        "../envs/distributions.yaml"
    script:
        "../../scripts/get_parameter_distributions.py"

init_dct = {str(k):v for k,v in enumerate(2**np.linspace(-0.4,3,50))}

rule grid_search:
    input:
        original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
        valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
    params:
        lambda wildcards: init_dct[wildcards.random_init] #'{random_init}' # delta umi mhc
    output:
        CAT_DIRECTORY + "/eval_clonotypes/grid_search/{random_init}.csv" #"grid_search.{umi_count_mhc_rel}.csv"
    script:
        "../../scripts/extract_concordance_table.py"
        
rule grid_search_swarm:
    input:
        CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
    params:
        '{random_init}' #?
    output:
        CAT_DIRECTORY + "/eval_clonotypes/grid_search/{random_init}.swarm.csv" #"grid_search.{umi_count_mhc_rel}.csv"
    script:
        "../../scripts/grid_search.swarm.py"
        
# #rule grid_search:
# #    input:
# #        original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
# #        valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
# #        min_thrs = CAT_DIRECTORY + "/eval_clonotypes/threshold_min.csv",
# #        var_dist = CAT_DIRECTORY + "/eval_clonotypes/param_dist.csv"
# #    params:
# #        random_init = '{random_init}',
# #        samples = 100000
# #    output:
# #        CAT_DIRECTORY + "/eval_clonotypes/grid_search.{random_init}.csv"
# #    script:
# #        "../../scripts/grid_search.py"

rule plot_grid:
    input:
        valid = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        grids = expand(CAT_DIRECTORY + "/eval_clonotypes/grid_search/{random_init}.csv", random_init=init_dct.keys()) #np.arange(100)
    output:
        plots = expand(PLT_DIRECTORY + "/eval_clonotypes/grid_search/opt_thr/grid.{ext}", ext=["pdf", "png"]),
        table = CAT_DIRECTORY + "/eval_clonotypes/threshold/opt.csv"
    script:
        "../../scripts/plot_grid.py"
        

 
#rule plot_thr_specificity_matrix:
#    input:
#        thresholds = expand(CAT_DIRECTORY + "/eval_clonotypes/threshold_{thrs_vals}.csv", thrs_vals=get_thresholds),
#        original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
#        valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
#    output:
#        plots = expand(PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{thrs_vals}/{{plot_type}}.{{data_type}}.{{thrs_type}}.{{ext}}", thrs_vals=get_thresholds)
#    params:
#        plot_type = lambda wildcards: wildcards.plot_type,
#        data_type = lambda wildcards: wildcards.data_type,
#        thrs_type = lambda wildcards: wildcards.thrs_type
#    conda:
#        "../envs/basic_dependencies.yaml"
#    priority:
#        90
#    script:
#        "../../scripts/E_plot_specificity_matrix.thr.py"
        
#rule plot_thr_specificity_matrix:
#    input:
#        thresholds = CAT_DIRECTORY + "/eval_clonotypes/threshold/{thrs_vals}.csv",
#        original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
#        valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
#    output:
#        plots = PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{thrs_vals}/{plot_type}.{data_type}.{thrs_type}.{ext}"
#    params:
#        plot_type = lambda wildcards: wildcards.plot_type,
#        data_type = lambda wildcards: wildcards.data_type,
#        thrs_type = lambda wildcards: wildcards.thrs_type
#    conda:
#        "../envs/basic_dependencies.yaml"
#    #priority:
#    #    90
#    script:
#        "../../scripts/E_plot_specificity_matrix.thr.py"

# Testing
rule plot_thr_specificity_matrix:
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
    #priority:
    #    90
    script:
        "../../scripts/E_plot_specificity_matrix-Copy1.thr.py"
        
rule plot_sample_specificity_matrix:
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
    input:
        CAT_DIRECTORY + "/similarity_assessment/{chain}_cdr3_all.txt"
    params:
        chain = CAT_DIRECTORY + "/similarity_assessment/{chain}"
    output:
        directory(CAT_DIRECTORY + "/similarity_assessment/chunks/{chain}/") # added chunks
    run:
        lst = np.loadtxt(input[0], dtype='str')
        # Write chunck files for parallel processing
        n = 200
        for i in range(0, len(lst), n):
            np.savetxt(output[0] + str(i) + '.txt', lst[i:i + n], fmt='%s') # permanently save
            
rule kernel_similarity_between_all_cdr3s:
    priority: 100
    input:
        s2s = config['s2s'], #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel",
        blf = config['blf'], #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50",
        qij = config['qij'], #"/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij",
        chu = CAT_DIRECTORY + "/similarity_assessment/chunks/{chain}/{chunk}.txt", # added chunks
        cdr = CAT_DIRECTORY + "/similarity_assessment/{chain}_cdr3_all.txt"
    params:
        chunk_i = '{chunk}',
        chain = '{chain}'
    output:
        chunk = CAT_DIRECTORY + "/similarity_assessment/chunks/{chain}/{chunk}.csv" # added chunks
    script:
        "../../scripts/K_kernel_similarity_between_all_cdr3s.py"
        
def aggregate_similarity_computations(wildcards):
    checkpoint_output = checkpoints.split_cdr3_seqs_in_chuncks.get(**wildcards).output[0]
    chunks = glob_wildcards(os.path.join(checkpoint_output, "{chunk}.txt")).chunk
    print('aggregate func')
    return expand(CAT_DIRECTORY + "/similarity_assessment/chunks/{chain}/{chunk}.csv", chain=wildcards.chain, chunk=chunks) # added chunks
        
rule aggregate_similarity_computations:
    input:
        aggregate_similarity_computations
    output:
        CAT_DIRECTORY + "/similarity_assessment/cdr3_{chain}.csv"
    shell:
        "cat {input} > {output}"
        
#TCR_DIRECTORY + '/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv' #CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"
#rule plot_tcr_kernel_similarity:
        
#rule plot_thr_similarity:
#    input:
#        threshold = CAT_DIRECTORY + "/eval_clonotypes/threshold/{thrs_vals}.csv", #opt
#        original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
#        valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
#        sim_tra = CAT_DIRECTORY + "/similarity_assessment/cdr3_a.csv",
#        sim_trb = CAT_DIRECTORY + "/similarity_assessment/cdr3_b.csv"
#    output:
#        plt_df = CAT_DIRECTORY + "/similarity_assessment/plt_df/{thrs_vals}.csv",
#        score_per_pep = PLT_DIRECTORY + "/eval_clonotypes/similarity/{thrs_vals}/score_per_pep.pdf",
#        score_pooled = PLT_DIRECTORY + "/eval_clonotypes/similarity/{thrs_vals}/score_pooled.pdf",
#        score_pooled_delta = PLT_DIRECTORY + "/eval_clonotypes/similarity/{thrs_vals}/score_pooled_delta.pdf",
#        auc = PLT_DIRECTORY + "/eval_clonotypes/similarity/{thrs_vals}/auc.pdf"
#    #priority:
#    #90
#    conda:
#        "../envs/basic_dependencies.yaml"
#    script:
#        "../../scripts/G_plot_tcr_kernel_similarity.py"

checkpoint index_similarity_per_peptide:
    input:
        valid = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        #sim_tra = CAT_DIRECTORY + "/similarity_assessment/cdr3_a.csv",
        #sim_trb = CAT_DIRECTORY + "/similarity_assessment/cdr3_b.csv",
        idx_tra = CAT_DIRECTORY + "/similarity_assessment/a_cdr3_all.txt",
        idx_trb = CAT_DIRECTORY + "/similarity_assessment/b_cdr3_all.txt"
    output:
        directory(CAT_DIRECTORY + "/similarity_assessment/tmp_cmd/") #split/map/
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/G_idx_sim_per_pep.py"
        
#rule split_similarity_per_peptide:
#    input:
#        CAT_DIRECTORY + "/similarity_assessment/tmp_cmd/{chain}.{pep}.sh"
#    output:
#        CAT_DIRECTORY + "/similarity_assessment/tmp_csv/{chain}.{pep}.csv"
#    shell:
#        """
#        sh {input} > {output}
#        """
        
rule split_similarity_per_peptide:
    #https://unix.stackexchange.com/questions/506207/fast-way-to-extract-lines-from-a-large-file-based-on-line-numbers-stored-in-anot/506226#506226
    #cut -f1 {input.idx} | while read rownum; do printf '%.12d\n' "$rownum"; done > {params.tmp}
    #join <(sort {params.tmp}) <(nl -w 12 -n rz {input.sim}) | cut -d ' ' -f 2- > {output}
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
    input:
        opt_thr = CAT_DIRECTORY + "/eval_clonotypes/threshold/opt.csv", #opt {thrs_vals}
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
    print(checkpoint_output)
    pep = glob_wildcards(os.path.join(checkpoint_output, "{chain}.{pep}.txt")).pep
    print(pep)
    print('aggregate func')
    return expand(CAT_DIRECTORY + "/similarity_assessment/sim_pep/{filtering_set}/{random_sample}/{pep}.csv",
                  filtering_set=wildcards.filtering_set, random_sample=wildcards.random_sample, pep=pep)
        
rule agg_similarity_per_peptide:
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
        expand(CAT_DIRECTORY + "/similarity_assessment/plt_df/{{filtering_set}}/{random_sample}.csv", random_sample=1) #np.arange(30)
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
    #priority:
    #90
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/G_plot_inter_v_inter_similarity.py"
        
rule get_filters:
    input:
        opt_thr = CAT_DIRECTORY + "/eval_clonotypes/threshold/opt.csv", #opt {thrs_vals}
        valid = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        gex = CAT_DIRECTORY + "/eval_clonotypes/threshold/gex.txt"
    params:
        flt = '{filtering_set}'
    output:
        flt = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.yaml",
        idx = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.csv"
    script:
        "../../scripts/H_set_filters.py"
        
rule filter_impact_staircase:
    input:
        df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
        flt = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.yaml",
        idx = CAT_DIRECTORY + "/eval_clonotypes/threshold/{filtering_set}.csv"
    output:
        touch(PLT_DIRECTORY + "/specificity_matrix/peptide_per_clonotype_by_gem_size/{filtering_set}/plots.done")
    conda:
        "../envs/basic_dependencies.yaml"
    script:
        "../../scripts/filtering_impact.staircase.py"
        
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
