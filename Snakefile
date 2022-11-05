workdir: "/home/tuba/herpov/tcr-pmhc-sc-project"
configfile: "snake_configs/config.yaml" #"config.yaml"
#include: "reports/Snakefile"

import numpy as np

MHC_DIRECTORY = config["parent_dir"] + "/data/" + config["exp"] + "_MHC_" + config["platform"]
TCR_DIRECTORY = config["parent_dir"] + "/data/" + config["exp"] + "_TCR"
CAT_DIRECTORY = config["parent_dir"] + "/data/" + config["exp"] + "_CAT_" + config["platform"] + "_" + config["mapping"] + "_" + config["barcode_system"]
PLT_DIRECTORY = config["parent_dir"] + "/plots/" + config["exp"] + "_" + config["platform"] + "_" + config["mapping"] + "_" + config["barcode_system"]

rule all:
	input:
		"initiation.done",#"test.done",
		#TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
		#MHC_DIRECTORY + "/reports/fastqc/" + config["platform"] + ".R1.gems.no_umi_fastqc.html",
		#TCR_DIRECTORY + "/reports/fastqc/TCR_S1_L001_R1_001_fastqc.html",
		expand("{DIRECTORY}/reports/gems/gem_counts.json", DIRECTORY=[MHC_DIRECTORY, TCR_DIRECTORY, CAT_DIRECTORY]),
		#MHC_DIRECTORY + "/reports/gems/gem_counts.longranger.lst",
		expand("{PLT_DIRECTORY}/reports/cutadapt_oligo_search/pie_{oligo}_present.pdf", PLT_DIRECTORY=PLT_DIRECTORY, oligo=['all', '1', '2', '3']),
		# PLT_DIRECTORY + "/reports/read_level/read_counts_mhc_pipeline.pdf",
		# PLT_DIRECTORY + "/reports/read_level/read_lengths_mhc_pipeline.pdf",
		# PLT_DIRECTORY + "/reports/read_level/umi_dist_tcr_pipeline.pdf",
		# PLT_DIRECTORY + "/reports/read_level/umi_dist_bc_pipeline.pdf",
		PLT_DIRECTORY + "/mapping/alignments.png",
		PLT_DIRECTORY + "/reports/venn/MHC_AB.pdf",
		# PLT_DIRECTORY + "/reports/annotation_dist/GEMs_per_barcode_template.pdf",
		PLT_DIRECTORY + "/reports/venn/TCR_MHC_AB.pdf",
		# PLT_DIRECTORY + "/reports/venn/TCR_BC_pipeline.pdf",
		PLT_DIRECTORY + "/reports/venn/HLA_concordance.pdf",
		expand("{DIRECTORY}/reports/confusion_multiplets/{min_bc}.pdf", DIRECTORY=PLT_DIRECTORY, min_bc=['b1','b2','b3','b4','b5']),
		# expand("{DIRECTORY}/reports/mhc_umi_dist_per_ct/ct_{top_cts}.pdf", DIRECTORY=PLT_DIRECTORY, top_cts=list(range(1, 11))),
		#CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
		#expand("{DIRECTORY}/specificity_matrix/{plot_type}/{ct_fmt}/{umi_delta}/{filtration}/{min_bc}.{min_tcr}.{ecs}.{ess}.pdf", DIRECTORY=PLT_DIRECTORY, plot_type='peptide_per_clonotype_by_gem_size', ct_fmt=['ct'], umi_delta='umi_delta0', filtration='no_filtration', min_bc='b1', min_tcr='t1', ecs='ecs_False', ess=['ess_False']),
		expand("{DIRECTORY}/specificity_matrix/{plot_type}/{ct_fmt}/{umi_delta}/{filtration}/{min_bc}.{min_tcr}.{ecs}.{ess}.pdf", DIRECTORY=PLT_DIRECTORY, plot_type=config['plot_type'], ct_fmt=['ct'], umi_delta=config['umi_delta'], filtration=config['filtration'], min_bc=config['bc_thresh'], min_tcr=config['tcr_thresh'], ecs='ecs_False', ess=['ess_False', 'ess_True']),
		expand("{DIRECTORY}/evaluation_per_filtration_level/{eval_type}.pdf", DIRECTORY=PLT_DIRECTORY, eval_type=['similarity','concordance']),
		expand("{DIRECTORY}/eval_clonotypes/specificity_matrix/{plot_type}.{data_type}.{thrs_type}.{ext}", DIRECTORY=PLT_DIRECTORY, plot_type=config['plot_type'], data_type=['train','test','total'], thrs_type=['raw','thr'], ext=['pdf','png'])

#if config["platform"] == "ILLUMINA":
#	rule merge_CD8_MHC_Illumina_reads:
#		input:
#			mhc = MHC_DIRECTORY + "/raw/" + config["raw_file"]["MHC"],
#			cd8 = MHC_DIRECTORY + "/raw/" + config["raw_file"]["CD8"]
#		output:
#			MHC_DIRECTORY + "/raw/" + config["platform"] + config["raw_file"]["EXTZ"]
#		run:
#			print(config["platform"], "Platform")
#			shell("./pipeline/barcodes/A0_merge_CD8_MHC.sh {input} {output}")

#rule rename_Illumina_headers:
#	input:
#		mhc = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/" + config["raw_file"]["MHC"],
#		cd8 = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/" + config["raw_file"]["CD8"]
#	output:
#		mhc = temp(config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/MHC-tmp" + config["raw_file"]["EXT"]),
#		cd8 = temp(config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/CD8-tmp" + config["raw_file"]["EXT"])
#	shell:
#		#""" zcat {input.mhc} | head | awk '{{print (NR%4 == 1) ? $0 "-MHC" : $0}}' """
#		""" zcat {input.mhc} | awk '{{print (NR%4 == 1) ? $0 "-MHC" : $0}}' > {output.mhc} """
#		""" zcat {input.cd8} | awk '{{print (NR%4 == 1) ? $0 "-CD8" : $0}}' | gzip -c > {output.cd8} """
#		#"mkdir -p raw"
#		#"mv {input} raw"
#
#rule merge_CD8_MHC_Illumina_reads:
#	input:
#		mhc = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/MHC-tmp" + config["raw_file"]["EXT"],
#		cd8 = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/CD8-tmp" + config["raw_file"]["EXT"]
#	output:
#		config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/" + config["platform"] + config["raw_file"]["EXT"]
#	shell:
#		"cat {input} > {output}"

rule print_experiment:
	output:
		temp(touch("initiation.done"))
	shell:
		"/home/tuba-nobackup/shared/ascii_art/bin/figlet {config[exp]}"

#import os
#if os.path.isfile(TCR_DIRECTORY + '/processed/mkfastq/index.csv'):
#	longranger_clean

rule make_mkfastq_index:
	output:
		TCR_DIRECTORY + '/processed/mkfastq_in/index.csv'
	shell:
		"echo Lane,Sample,Index > {output}; echo 1,TCR,{config[index]} >> {output}"

rule run_cellranger_mkfastq:
	input:
		TCR_DIRECTORY + '/raw',
		TCR_DIRECTORY + '/processed/mkfastq_in/index.csv'
	output:
		TCR_DIRECTORY + '/processed/mkfastq_out/Stats/DemuxSummaryF1L1.txt'
	shell:
		'./pipeline/tcr/A_run_cellranger_mkfastq.sh {input} {output}'

rule run_cellranger_vdj:
	input:
		#TCR_DIRECTORY + '/processed/mkfastq_out'
		TCR_DIRECTORY + '/processed/mkfastq_out/Stats/DemuxSummaryF1L1.txt'
	output:
		TCR_DIRECTORY + '/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv'
	shell:
		'./pipeline/tcr/B_run_cellranger_vdj.sh {input} {output}'

rule extract_GEMs:
	input:
		TCR_DIRECTORY + '/processed/mkfastq_out/Stats/DemuxSummaryF1L1.txt'
	output:
		TCR_DIRECTORY + "/reports/gems/gem_counts.raw.lst",
		TCR_DIRECTORY + "/reports/umis/umi_per_gem.raw.tab"
	shell:
		"./tools/extract_reads_per_GEM_TCR.sh {input} {output}"

rule clean_augment_tcr:
	input:
		TCR_DIRECTORY + '/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv',
		TCR_DIRECTORY + '/processed/cellranger_out/TCR_VDJ/outs/consensus_annotations.csv'
	output:
		TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
		TCR_DIRECTORY + "/reports/gems/gem_counts.json"
	script:
		"./scripts/00_clean_augment_tcr.py"

if config["platform"] == "IONTORRENT":
	rule A1_cut_adapter:
		input:
			directory = MHC_DIRECTORY + "/raw/IONTORRENT_S1_L001_R1_001.fastq" # added filename!
		output:
			filename = MHC_DIRECTORY + "/processed/cut_adapter/" + config["platform"] + "_" + config["raw_file"]["EXT"]
		shell:
			"./pipeline/barcodes/A1_cut_adapter.sh {input.directory} {output.filename} {config[platform]}"

	rule A2_run_longranger:
		params:
			out_dir = MHC_DIRECTORY + "/processed/longranger_out"
		input:
			filename = MHC_DIRECTORY + "/processed/cut_adapter/" + config["platform"] + "_" + config["raw_file"]["EXT"]
		output:
			filename = MHC_DIRECTORY + "/processed/longranger_out/" + config["platform"] + "/outs/barcoded.fastq.gz"
			#directory = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/processed/longranger_out"
		shell:
			"./pipeline/barcodes/A2_run_longranger.sh {input.filename} {output.filename} {params.out_dir} {config[platform]}"

	rule B_clean_longranger_output:
		params:
			in_dir = MHC_DIRECTORY + "/processed/longranger_out",
			out_dir = MHC_DIRECTORY + "/processed/longranger_clean"
		input:
			filename = MHC_DIRECTORY + "/processed/longranger_out/" + config["platform"] + "/outs/barcoded.fastq.gz"
		output:
			filename = MHC_DIRECTORY + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fq"
		shell:
			"./pipeline/barcodes/B_clean_longranger_output.sh {input.filename} {params.in_dir} {output.filename} {params.out_dir} {config[platform]}"

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

# OBS! Update all oligo_a.fa files to also contain the A4000 sequence!
if config["barcode_system"] == "AKB":
	rule A_construct_template_database:
		input:
			oligo_a = MHC_DIRECTORY + "/barcode_library/oligo_a.fa",
			oligo_b = MHC_DIRECTORY + "/barcode_library/oligo_b.fa",
			samples = MHC_DIRECTORY + "/barcode_library/sample.fa",
			nomencl = MHC_DIRECTORY + "/barcode_library/barcode_specificity_annotations.xlsx"
		output:
			MHC_DIRECTORY + "/barcode_library/barcode_templates.fa"
		conda:
			"snake_envs/biopython.yaml"
		script:
			"./scripts/A0_construct_template_database.py"
else:
	rule B_construct_template_database:
		input:
			oligo_1os = MHC_DIRECTORY + "/barcode_library/oligo_1os.fa",
			oligo_mhc = MHC_DIRECTORY + "/barcode_library/oligo_mhc.fa",
			oligo_cdx = MHC_DIRECTORY + "/barcode_library/oligo_cdx.fa",
			nomencl = MHC_DIRECTORY + "/barcode_library/barcode_specificity_annotations.xlsx"
		output:
			MHC_DIRECTORY + "/barcode_library/barcode_templates.fa"
		conda:
			"snake_envs/biopython.yaml"
		script:
			"./scripts/A1_construct_template_database.py"

# Replaced by the two following rules
# rule E_kma_mapping:
# 	params:
# 		arguments = config["mapping_args"],
# 		#out_dir = MHC_DIRECTORY + "/kma/kma" + config["kma-args"].replace(" ", "") # arguments are added to the dirname in the shell file
# 	input:
# 		templates = MHC_DIRECTORY + "/barcode_library/barcode_templates.fa",
# 		queries = MHC_DIRECTORY + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fq"
# 	output:
# 		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.frag.gz"
# 	shell:
# 		"./pipeline/barcodes/E_kma_mapping.sh {input.templates} {input.queries} {params.arguments} {output}"

rule kma_index:
	input:
		MHC_DIRECTORY + "/barcode_library/barcode_templates.fa"
	output:
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/templates.name"
	params:
		prefix = lambda wildcards, output: os.path.splitext(output[0])[0]
	shell:
		"/home/tuba-nobackup/shared/src/kma/kma index -i {input} -o {params.prefix}"

rule kma_map:
	input:
		reads = MHC_DIRECTORY + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fq",
		templates = MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/templates.name"
	output:
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.frag.gz"
	params:
		input_prefix = lambda wildcards, input: os.path.splitext(input[1])[0],
		output_prefix = lambda wildcards, output: output[0].rsplit('.', 2)[0],
		arguments = config["mapping_args"]
	shell:
		"/home/tuba-nobackup/shared/src/kma/kma -i {input.reads} -o {params.output_prefix} -t_db {params.input_prefix} {params.arguments}"

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

rule plot_venn_MHC_AB:
	input:
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean." + config['barcode_system'] + ".augmented.gz"
	output:
		PLT_DIRECTORY + "/reports/venn/MHC_AB.pdf"
	conda:
		"snake_envs/venn.yaml"
	script:
		"./scripts/plot_venn_MHC_AB.py"

if config["mapping"] == "KMA":
	rule clean_kma:
		input:
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.frag.gz"
		output:
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean.gz"
		script:
			"./scripts/A_clean_kma.py"

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

	# --cores 20
	rule extract_umis:
		input:
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.frag.gz",
			MHC_DIRECTORY + "/barcode_library/sample.fa",
			MHC_DIRECTORY + "/barcode_library/oligo_a.fa",
			MHC_DIRECTORY + "/barcode_library/oligo_b.fa",
			MHC_DIRECTORY + "/barcode_library/barcode-information.fa"
		output:
			MHC_DIRECTORY + "/barcode_library/umi.tsv"
		#conda:
		#	"snake_envs/umi.yaml"
		shell:
			"/home/tuba-nobackup/shared/R/R-3.6.1/bin/Rscript /home/tuba/kamilla/10x-barcoding/scripts/parse-kma-results.R {input} {output}"

	rule augment_mapping:
		params:
			config['barcode_system']
		input:
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean.gz",
			MHC_DIRECTORY + "/barcode_library/barcode_specificity_annotations.xlsx",
			MHC_DIRECTORY + "/barcode_library/detected_responses_annotation.xlsx",
			MHC_DIRECTORY + "/barcode_library/umi.tsv"
		output:
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean." + config['barcode_system'] + ".augmented.gz",
			MHC_DIRECTORY + "/reports/gems/gem_counts.json"
		script:
			"./scripts/B_augment_mapping.py"

	rule comb_barcodes_TCR:
		params:
			config['barcode_system']
		input:
			TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean." + config['barcode_system'] + ".augmented.gz"
		output:
			CAT_DIRECTORY + "/tables/tcr_barcode.csv"
		script:
			"./scripts/C_comb_barcodes_TCR.py"

elif config["mapping"] == "BLAST":
	rule augment_mapping:
		params:
			config['barcode_system']
		input:
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + "/blast.annotated.clean.tsv",
			MHC_DIRECTORY + "barcode_library/barcode_specificity_annotations.tab",
			MHC_DIRECTORY + "barcode_library/detected_responses_annotation.xlsx"
		output:
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + "/blast.annotated.clean." + config['barcode_system'] + ".augmented.tsv"
		script:
			"./scripts/B_augment_mapping.py"

	rule comb_barcodes_TCR:
		params:
			config['barcode_system']
		input:
			TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv",
			MHC_DIRECTORY + "/mapping/" + config["mapping"] + "/blast.annotated.clean." + config['barcode_system'] + ".augmented.tsv"
		output:
			CAT_DIRECTORY + "/tables/tcr_barcode.csv"
		script:
			"./scripts/C_comb_barcodes_TCR.py"

rule clean_tcr_barcodes: # RENAME: clean_tcr_barcodes
	input:
		CAT_DIRECTORY + "/tables/tcr_barcode.csv",
		"/tools/VDJdb.csv"
	params:
		TCR_DIRECTORY + "/library/clone_rearrangement.tsv"
	output:
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
		CAT_DIRECTORY + "/reports/gems/gem_counts.json"
	script:
		"./scripts/D_clean_tcr_barcodes.py"

#rule imputations: # Not yet tested
#	input:
#		df = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"
#	output:
#		output = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.imputed.csv",
#		report = CAT_DIRECTORY + "/reports/imputations.csv"
#	script:
#		"./scripts/F_imputations.py"

rule eval_clonotypes: # Not yet tested
	input:
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"
	params:
		PLT_DIRECTORY + "/eval_clonotypes/%s/%d.pdf"
	output:
		CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv",
		touch(expand(PLT_DIRECTORY + "/eval_clonotypes/{flag}/dir.done", flag=['significant','insignificant']))
	script:
		"./scripts/F_comp_cred_specificities.py"

rule grid_search:
	input:
		original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
		valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
	output:
		expand(CAT_DIRECTORY + "/eval_clonotypes/grid_search.{umi_count_mhc_rel}.csv", umi_count_mhc_rel=list((np.linspace(0,0.2,20) * 10000).astype(int).astype(str)))
	script:
		"./scripts/extract_concordance_table.py"

rule plot_grid:
	input:
		CAT_DIRECTORY + "/eval_clonotypes/grid_search.0.csv"
	output:
		PLT_DIRECTORY + "/eval_clonotypes/opt_thr/grid.pdf",
		PLT_DIRECTORY + "/eval_clonotypes/opt_thr/threshold.csv"
	script:
		"./scripts/plot_grid.py"

rule plot_thr_specificity_matrix:
	input:
		threshold = PLT_DIRECTORY + "/eval_clonotypes/opt_thr/threshold.csv",
		original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
		valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
	output:
		PLT_DIRECTORY + "/eval_clonotypes/specificity_matrix/{plot_type}.{data_type}.{thrs_type}.{ext}"
	priority:
		90
	script:
		"./scripts/E_plot_specificity_matrix.thr.py"

rule plot_thr_similarity:
	input:
		threshold = PLT_DIRECTORY + "/eval_clonotypes/opt_thr/threshold.csv",
		original = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
		valid_df = CAT_DIRECTORY + "/eval_clonotypes/valid_ct.csv"
	output:
		score_per_pep = PLT_DIRECTORY + "/eval_clonotypes/similarity/score_per_pep.pdf",
		score_pooled = PLT_DIRECTORY + "/eval_clonotypes/similarity/score_pooled.pdf",
		score_pooled_delta = PLT_DIRECTORY + "/eval_clonotypes/similarity/score_pooled_delta.pdf"
	priority:
		90
	script:
		"./scripts/G_plot_similarity.py"

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

rule plot_venn_TCR_MHC_AB:
	input:
		MHC_DIRECTORY + "/mapping/" + config["mapping"] + config["mapping_args"].replace(" ", "") + "/output/mapping.clean." + config['barcode_system'] + ".augmented.gz",
		TCR_DIRECTORY + "/augmented/tcr.clean.augmented.csv"
	output:
		PLT_DIRECTORY + "/reports/venn/TCR_MHC_AB.pdf"
	conda:
		"snake_envs/venn.yaml"
	priority:
	    100
	script:
		"./scripts/plot_venn_TCR_MHC_AB.py"

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

rule plot_venn_HLA_concordance:
	input:
		CAT_DIRECTORY + "/tables/tcr_barcode.csv",
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.imputed.csv"
	output:
		PLT_DIRECTORY + "/reports/venn/HLA_concordance.pdf"
	conda:
		"snake_envs/venn.yaml"
	script:
		"./scripts/plot_venn_HLA_concordance.py"

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


rule plot_confusion_multiplets:
	input:
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv" # Not imputed?
	output:
		PLT_DIRECTORY + "/reports/confusion_multiplets/{min_bc}.pdf"
	script:
		"./scripts/plot_confusion_multiplets.py"

#print(expand("{DIRECTORY}/specificity_matrix/{plot_type}/{{ct_fmt}}/{{umi_delta}}/{{filtration}}/{{min_bc}}.{{min_tcr}}.{{ecs}}.{{ess}}.pdf", DIRECTORY=PLT_DIRECTORY, plot_type=config['plot_type']))
rule plot_specificity_matrix:
	params:
		config['version']
	input:
		CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv" #.imputed
	output:
		#plot = directory(expand("{DIRECTORY}/specificity_matrix/{plt}/{ct_fmt}/{filtration}/",
		#						DIRECTORY=PLT_DIRECTORY, plt='peptide_per_clonotype_by_gem_size', ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs'])),
		#report = directory(expand("{DIRECTORY}/specificity_matrix/{plt}/{ct_fmt}/{filtration}/{type}/",
		#						DIRECTORY=CAT_DIRECTORY, plt='peptide_per_clonotype_by_gem_size', ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs'], type=['unique_tcrs', 'unique_gems', 'concordance_comparison_matrix']))
		# Other option..
		#plot = PLT_DIRECTORY + "/specificity_matrix/{plot_type}/{ct_fmt}/{umi_delta}/{filtration}/{min_bc}.{min_tcr}.{ecs}.{ess}.pdf", #
		plots = expand("{DIRECTORY}/specificity_matrix/{plot_type}/{{ct_fmt}}/{{umi_delta}}/{{filtration}}/{{min_bc}}.{{min_tcr}}.{{ecs}}.{{ess}}.pdf", DIRECTORY=PLT_DIRECTORY, plot_type=config['plot_type']),
		#expand("{DIRECTORY}/specificity_matrix/{plt}/{ct_fmt}/{filtration}/{min_bc}.{min_tcr}.{ecs}.{ess}.pdf", DIRECTORY=PLT_DIRECTORY, plt='peptide_per_clonotype_by_gem_size', ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs'], min_bc=['b1','b2','b3','b4','b5'], min_tcr=['t1','t2','t4','t6','t8','t10'], ecs='ecs_False', ess='ess_False')
		report = expand("{DIRECTORY}/specificity_matrix/{{ct_fmt}}/{{umi_delta}}/{{filtration}}/{unique_list}/{{min_bc}}.{{min_tcr}}.{{ecs}}.{{ess}}.lst", DIRECTORY=CAT_DIRECTORY, unique_list=['unique_tcrs', 'unique_gems']), #peptide_per_clonotype_by_gem_size/
		# expand("{DIRECTORY}/specificity_matrix/peptide_per_clonotype_by_gem_size/{ct_fmt}/{filtration}/{unique_list}/{min_bc}.{min_tcr}.{ecs}.{ess}.lst", DIRECTORY=CAT_DIRECTORY, ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs'], unique_list=['unique_tcrs', 'unique_gems'], min_bc=['b1','b2','b3','b4','b5'], min_tcr=['t1','t2','t4','t6','t8','t10'], ecs='ecs_False', ess='ess_False'),
		#report_concordance = expand("{DIRECTORY}/specificity_matrix/peptide_per_clonotype_by_gem_size/{{ct_fmt}}/{{filtration}}/concordance_comparison_matrix/{{umi_filtrations}}/{coordinate_file}.lst", DIRECTORY=CAT_DIRECTORY, coordinate_file=['X', 'Y', 'Z_conc', 'Z_gems'])
		# expand("{DIRECTORY}/specificity_matrix/peptide_per_clonotype_by_gem_size/{ct_fmt}/{filtration}/concordance_comparison_matrix/{coordinate_file}.lst", DIRECTORY=CAT_DIRECTORY, ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs'], coordinate_file=['X', 'Y', 'Z_conc', 'Z_gems'])

	conda:
		"snake_configs/env_ipynb.yaml"
	priority:
		90
	script:
		"./scripts/E_plot_specificity_matrix.py"

# https://groups.google.com/forum/#!topic/snakemake/bYvHjE7YA5I
# https://bitbucket.org/snakemake/snakemake/issues/790/explaining-expand-and-variable-to-mask
CT_FMT, UMI_DELTA, FILTRATION, UNIQUE_GEM_FILENAMES = glob_wildcards(CAT_DIRECTORY + "/specificity_matrix/{ct_fmt}/{umi_delta}/{filtration}/unique_gems/{filename}.ecs_False.ess_False.lst")
CT_FMT, UMI_DELTA, FILTRATION = list(set(CT_FMT)), list(set(UMI_DELTA)), list(set(FILTRATION))
# OBS! I am only looking at false false cases!
#print(UNIQUE_GEM_FILENAMES)

rule kernel_similarity_between_all_cdr3s:
	input:
		s2s = "/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel",
		blf = "/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50",
		qij = "/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij",
		df = TCR_DIRECTORY + '/processed/cellranger_out/TCR_VDJ/outs/all_contig_annotations.csv' #CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv"
	output:
		A = CAT_DIRECTORY + "/similarity_assessment/cdr3/kernel_similarity_scores/cdr3_a.csv",
		B = CAT_DIRECTORY + "/similarity_assessment/cdr3/kernel_similarity_scores/cdr3_b.csv"
	script:
		"./scripts/kernel_similarity_between_all_cdr3s.py"

rule plot_tcr_kernel_similarity:
	input:
		#df = CAT_DIRECTORY + "/similarity_assessment/cdr3/kernel_similarity_scores/cdr3.csv",
		df = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv", #.imputed
		sim_tra = CAT_DIRECTORY + "/similarity_assessment/cdr3/kernel_similarity_scores/cdr3_a.csv",
		sim_trb = CAT_DIRECTORY + "/similarity_assessment/cdr3/kernel_similarity_scores/cdr3_b.csv",
		unique_gems = expand("{DIRECTORY}/specificity_matrix/{{ct_fmt}}/{{umi_delta}}/{{filtration}}/unique_gems/{umi_filtrations}.ecs_False.ess_False.lst",
							  DIRECTORY=CAT_DIRECTORY, umi_filtrations=UNIQUE_GEM_FILENAMES), #ct_fmt=CT_FMT, filtration=FILTRATION,
		peptides = MHC_DIRECTORY + "/barcode_library/barcode_specificity_annotations.xlsx"
		#unique_gems = expand("{DIRECTORY}/specificity_matrix/peptide_per_clonotype_by_gem_size/{{ct_fmt}}/{{filtration}}/unique_gems/{{umi_filtrations}}.lst",
		#					  DIRECTORY=CAT_DIRECTORY, umi_filtrations=UNIQUE_GEM_FILENAMES) 
	output:
		#plot = directory(expand("{DIRECTORY}/similarity_assessment/cdr3/{{ct_fmt}}/{{filtration}}/{plt}/{min_bc}/{min_tcr}/",
		#						DIRECTORY=PLT_DIRECTORY,
		#						plt=['pooled'], #'individual',
		#						min_bc=['b1', 'b2', 'b3', 'b4', 'b5'],
		#						min_tcr=['t1', 't2', 't4', 't6', 't8', 't10'])), #ct_fmt=CT_FMT,filtration=FILTRATION,
		plot   = expand("{DIRECTORY}/similarity_assessment/cdr3/{{ct_fmt}}/{{umi_delta}}/{{filtration}}/pooled/{plot_type}/{umi_filtrations}.pdf",
						  DIRECTORY=PLT_DIRECTORY, plot_type=['boxplot', 'pieplot'], umi_filtrations=UNIQUE_GEM_FILENAMES),
		file_X = CAT_DIRECTORY + "/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/X",
		file_Y = CAT_DIRECTORY + "/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/Y",
		siml_Z = CAT_DIRECTORY + "/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/Z_siml",
		gems_Z = CAT_DIRECTORY + "/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/Z_gems"

	script:
		"./scripts/plot_tcr_kernel_similarity.py"

rule plot_similarity_per_filtration_level:
	input:
		table  = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv", #.imputed
		file_X = expand("{DIRECTORY}/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/X",
					DIRECTORY=CAT_DIRECTORY, ct_fmt=CT_FMT, umi_delta=UMI_DELTA, filtration=FILTRATION), #['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs']),
		file_Y = expand("{DIRECTORY}/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/Y",
					DIRECTORY=CAT_DIRECTORY, ct_fmt=CT_FMT, umi_delta=UMI_DELTA, filtration=FILTRATION), #ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs']),
		siml_Z = expand("{DIRECTORY}/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/Z_siml",
					DIRECTORY=CAT_DIRECTORY, ct_fmt=CT_FMT, umi_delta=UMI_DELTA, filtration=FILTRATION), #ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs']),
		gems_Z = expand("{DIRECTORY}/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/Z_gems",
					DIRECTORY=CAT_DIRECTORY, ct_fmt=CT_FMT, umi_delta=UMI_DELTA, filtration=FILTRATION) #ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs'])
	output:
		PLT_DIRECTORY + "/evaluation_per_filtration_level/similarity.pdf"
	script:
		"./scripts/plot_similarity_per_filtration_level.py"

# Rule not of interest
rule plot_concordance_per_filtration_level:
	input:
		table = CAT_DIRECTORY + "/tables/tcr_barcode.cleaned.csv",
		file_X = expand("{DIRECTORY}/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/X",
					DIRECTORY=CAT_DIRECTORY, ct_fmt=CT_FMT, umi_delta=UMI_DELTA, filtration=FILTRATION), #['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs']),
		file_Y = expand("{DIRECTORY}/similarity_assessment/cdr3/similarity_comparison_matrix/{ct_fmt}/{umi_delta}/{filtration}/Y",
					DIRECTORY=CAT_DIRECTORY, ct_fmt=CT_FMT, umi_delta=UMI_DELTA, filtration=FILTRATION), #ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs']),
		#conc_Z = expand("{DIRECTORY}/specificity_matrix/peptide_per_clonotype_by_gem_size/{ct_fmt}/{filtration}/concordance_comparison_matrix/Z_conc",
		#			DIRECTORY=CAT_DIRECTORY, ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs']),
		#gems_Z = expand("{DIRECTORY}/specificity_matrix/peptide_per_clonotype_by_gem_size/{ct_fmt}/{filtration}/concordance_comparison_matrix/Z_gems",
		#			DIRECTORY=CAT_DIRECTORY, ct_fmt=['ct', 'num_clonotype'], filtration=['no_filtration', 'exclude_single-chain_TCRs', 'exclude_ambiguous_and_single-chain_TCRs'])
	output:
		PLT_DIRECTORY + "/evaluation_per_filtration_level/concordance.pdf"
	script:
		"./scripts/plot_concordance_per_filtration_level.py"


