rule C_run_fastqc_mhc:
	input:
		MHC_DIRECTORY + "/raw/" + config["platform"] + "_" + config["raw_file"]["EXTZ"],
		MHC_DIRECTORY + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fq"
	output:
		MHC_DIRECTORY + "/reports/fastqc/" + config["platform"] + "_S1_L001_R1_001_fastqc.html"
	shell:
		"./pipeline/barcodes/D_run_fastqc.sh {output} {input}"

rundependent_dir, = glob_wildcards(TCR_DIRECTORY + "/processed/mkfastq_out/{dir}/TCR")
#rule C_run_fastqc_tcr:
#	input:
#		TCR_DIRECTORY + "/processed/mkfastq_out/%s/TCR/TCR_S1_L001_I1_001.fastq.gz" %rundependent_dir[0],
#		TCR_DIRECTORY + "/processed/mkfastq_out/%s/TCR/TCR_S1_L001_R1_001.fastq.gz" %rundependent_dir[0],
#		TCR_DIRECTORY + "/processed/mkfastq_out/%s/TCR/TCR_S1_L001_R2_001.fastq.gz" %rundependent_dir[0]
#	output:
#		TCR_DIRECTORY + "/reports/fastqc/TCR_S1_L001_R1_001_fastqc.html"
#	shell:
#		"./pipeline/barcodes/D_run_fastqc.sh {output} {input}"

################################################################################
rule extract_GEMs:
	input:
		TCR_DIRECTORY + '/processed/mkfastq_out/Stats/DemuxSummaryF1L1.txt'
	output:
		TCR_DIRECTORY + "/reports/gems/gem_counts.raw.lst",
		TCR_DIRECTORY + "/reports/umis/umi_per_gem.raw.tab"
	shell:
		"./tools/extract_reads_per_GEM_TCR.sh {input} {output}"

rule count_gems_from_raw: # Perhaps change this to extract putative GEMs after adapter removal.
	input:
		MHC_DIRECTORY + "/raw/" + config["platform"] + "_" + config["raw_file"]["EXT"] 
	output:
		MHC_DIRECTORY + "/reports/gems/gem_counts.raw.lst"
	params:
		adapter_5p = "CACGACGCTCTTCCGATCT"
	shell:
		"""grep -E '^{params.adapter_5p}[ATCG]{{20}}' {input} | cut -c20-35 | awk '{{print $0"-1"}}' | sort | uniq > {output}"""

rule count_umis_from_raw: # Perhaps change this to extract putative GEMs after adapter removal.
	input:
		MHC_DIRECTORY + "/raw/" + config["platform"] + "_" + config["raw_file"]["EXT"] 
	output:
		MHC_DIRECTORY + "/reports/umis/umi_per_gem.raw.tab"
	params:
		adapter_5p = "CACGACGCTCTTCCGATCT"
	shell:
		"""grep -E '^{params.adapter_5p}[ATCG]{{20}}' {input} | cut -c20-45 | sort | uniq | cut -c20-35 | awk '{{print $0"-1"}}' | sort | uniq -c | awk '{{print $2, $1}}' > {output}"""

rule count_gems_from_longranger:
	input:
		MHC_DIRECTORY + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fq"
	output:
		MHC_DIRECTORY + "/reports/gems/gem_counts.longranger.lst"
	shell:
		"grep 'BX:Z:' {input} | cut -d':' -f5 | sort | uniq > {output}"

###########################################################################
rule extract_nonoverlapping_TCR_GEMs:
	input:
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.frag.gz",
		TCR_DIRECTORY + '/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv',
		TCR_DIRECTORY + '/processed/mkfastq_out/Stats/DemuxSummaryF1L1.txt'
	output:
		TCR_DIRECTORY + "/reports/umis/umi_per_gem.filtered.tab",
		TCR_DIRECTORY + "/reports/gems/gem_counts.filtered.lst"
	shell:
		"./tools/extract_non-overlapping_TCR_GEMs.sh {input} {output}"

# MHC Reports
rule count_reads:
	input:
		MHC_DIRECTORY + "/raw/" + config["platform"] + "_" + config["raw_file"]["EXT"], #EXTZ
		MHC_DIRECTORY + "/processed/cut_adapter/" + config["platform"] + "_" + config["raw_file"]["EXT"],
		MHC_DIRECTORY + "/processed/longranger_out/" + config["platform"] + "/outs/barcoded.fastq.gz",
		MHC_DIRECTORY + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fq",
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.frag.gz"
	output:
		MHC_DIRECTORY + "/reports/read_level/read_counts_mhc_pipeline.txt"
	shell:
		"./tools/count_reads.sh {input} {output}"

rule plot_read_counts_mhc_pipeline:
	input:
		MHC_DIRECTORY + "/reports/read_level/read_counts_mhc_pipeline.txt"
	output:
		PLT_DIRECTORY + "/reports/read_level/read_counts_mhc_pipeline.pdf"
	script:
		"./scripts/plot_read_counts_mhc_pipeline.py"

rule count_read_lengths:
	params:
		MHC_DIRECTORY + "/reports/read_level/read_lengths"
	input:
		MHC_DIRECTORY + "/raw/" + config["platform"] + "_" + config["raw_file"]["EXT"], #Z
		MHC_DIRECTORY + "/processed/cut_adapter/" + config["platform"] + "_" + config["raw_file"]["EXT"],
		MHC_DIRECTORY + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fq",
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.frag.gz"
	output:
		MHC_DIRECTORY + "/reports/read_level/read_lengths/raw.txt",
		MHC_DIRECTORY + "/reports/read_level/read_lengths/cut_adapter.txt",
		MHC_DIRECTORY + "/reports/read_level/read_lengths/longranger_clean.txt",
		MHC_DIRECTORY + "/reports/read_level/read_lengths/output.txt"
	shell:
		"./tools/count_read_lengths.sh {input} {output} {params}"

rule plot_read_lengths_mhc_pipeline:
	input:
		MHC_DIRECTORY + "/reports/read_level/read_lengths/raw.txt",
		MHC_DIRECTORY + "/reports/read_level/read_lengths/cut_adapter.txt",
		MHC_DIRECTORY + "/reports/read_level/read_lengths/longranger_clean.txt",
		MHC_DIRECTORY + "/reports/read_level/read_lengths/output.txt"
	output:
		PLT_DIRECTORY + "/reports/read_level/read_lengths_mhc_pipeline.pdf"
	script:
		"./scripts/plot_read_lengths_mhc_pipeline.py"

rule cutadapt_oligo_search:
	input:
		MHC_DIRECTORY + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fq"
	output:
		MHC_DIRECTORY + "/reports/cutadapt_oligo_search/headers/oligo_annotated_headers.lst",
		MHC_DIRECTORY + "/reports/cutadapt_oligo_search/headers/oligo_annotated_headers.uniq_count.lst"
	shell:
		"./tools/cutadapt_oligo_search.sh {input} {output}"

rule plot_cutadapt_oligo_search:
	params:
		PLT_DIRECTORY + "/reports/cutadapt_oligo_search/"
	input:
		MHC_DIRECTORY + "/reports/cutadapt_oligo_search/headers/oligo_annotated_headers.uniq_count.lst"
	output:
		expand("{PLT_DIRECTORY}/reports/cutadapt_oligo_search/pie_{oligo}_present.pdf", PLT_DIRECTORY=PLT_DIRECTORY, oligo=['all', '1', '2', '3'])
	script:
		"./scripts/plot_cutadapt_oligo_search.py"

################################################################
rule plot_kma_mapping_stats:
	input:
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean.gz"
	output:
		PLT_DIRECTORY + "/mapping/uncertainty_dist.pdf",
		PLT_DIRECTORY + "/mapping/reads_per_template.pdf",
		PLT_DIRECTORY + "/mapping/gems_per_template.pdf",
		PLT_DIRECTORY + "/mapping/reads_per_gem.pdf",
		PLT_DIRECTORY + "/mapping/alignments.png"
	script:
		"./scripts/plot_kma_mapping_stats.py"

#################################################################
rule plot_annotations_dist:
	input:
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.frag.gz",
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean." + config['barcode_system'] + ".augmented.gz",
		TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
		MHC_DIRECTORY + "/barcode_library/barcode_specificity_annotations.xlsx"
	output:
		PLT_DIRECTORY + "/reports/annotation_dist/GEMs_per_barcode_template.pdf",
		PLT_DIRECTORY + "/reports/annotation_dist/GEMs_per_clonotype.pdf"
	script:
		"./scripts/plot_annotation_dist.py"

###################################################################
rule plot_venn_GEMs_throughout_pipeline:
	input:
		BAR_GEMS = expand("{DIRECTORY}/reports/gems/gem_counts.{suffix}", DIRECTORY=MHC_DIRECTORY, suffix=['raw.lst','longranger.lst','json']),
		TCR_GEMS = expand("{DIRECTORY}/reports/gems/gem_counts.{suffix}", DIRECTORY=TCR_DIRECTORY, suffix=['raw.lst','filtered.lst','json']) #TCR_DIRECTORY + "/reports/gems/gem_counts.json"
	output:
		PLT_DIRECTORY + "/reports/venn/TCR_BC_pipeline.pdf"
	conda:
		"snake_envs/venn.yaml"
	script:
		"./scripts/plot_venn_GEMs_throughout_pipeline.py"

##################################################################
rule plot_TCR_UMI_distribution:
	input:
		TCR_DIRECTORY + "/reports/umis/umi_per_gem.raw.tab",
		TCR_DIRECTORY + "/reports/umis/umi_per_gem.filtered.tab",
		TCR_DIRECTORY + '/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv',
		TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"
	output:
		PLT_DIRECTORY + "/reports/read_level/umi_dist_tcr_pipeline.pdf"
	script:
		"./scripts/plot_tcr_umi_distribution.py"

rule plot_BC_UMI_distribution:
	input:
		MHC_DIRECTORY + "/reports/umis/umi_per_gem.raw.tab",
		MHC_DIRECTORY + "/barcode_library/umi.tsv",
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean." + config['barcode_system'] + ".augmented.gz",
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"
	output:
		PLT_DIRECTORY + "/reports/read_level/umi_dist_bc_pipeline.pdf"
	script:
		"./scripts/plot_bc_umi_distribution.py"

rule plot_MHC_UMI_dist_per_clonotype:
	input:
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv" #.imputed
	output:
		PLT_DIRECTORY + "/reports/mhc_umi_dist_per_ct/ct_{top_cts}.pdf"
	script:
		"./scripts/plot_mhc_umi_dist_per_clonotype.py"

####################################################################

