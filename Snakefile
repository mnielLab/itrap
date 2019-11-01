# workdir: /home/tuba/herpov/tcr-pmhc-sc-project
configfile: "config.yaml"

rule all:
	input:
		config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fa"
		#expand("{PARENT_DIR}/data/{EXP}_{PLATFORM}/raw/{PLATFORM}{FILE_EXT}",
		#	PARENT_DIR=config["parent_dir"], EXP=config["exp"], PLATFORM=config["platform"], FILE_EXT=config["raw_file_ext"])
		#"/home/tuba/herpov/tcr-pmhc-sc-project/data/exp2_MHC_ILLUMINA/raw/ILLUMINA-2_S1_L001_R1_001.fastq.gz"

rule merge_CD8_MHC_Illumina_reads:
	input:
		mhc = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/" + config["raw_file"]["MHC"],
		cd8 = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/" + config["raw_file"]["CD8"]
	output:
		config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw/" + config["platform"] + config["raw_file"]["EXTZ"]
	shell:
		"./pipeline/barcodes/A0_merge_CD8_MHC.sh {input} {output}"
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

rule cutadapt_reads:
	input:
		directory = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/raw"
	output:
		directory = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/processed/cut_adapter"
	shell:
		".pipeline/barcodes/A1_cut_adapter.sh {input.directory} {output.directory} {config[platform]} {config[raw_file][EXT]}"

rule longranger_annotates_gems:
	input:
		filename = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/processed/cut_adapter/" + config["platform"] + "_" + config["raw_file"]["EXT"]
	output:
		directory = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/processed/longranger_out"
	shell:
		".pipeline/barcodes/A2_run_longranger.sh {input.filename} {output.directory} {config[platform]}"

rule clean_longranger_output:
	input:
		directory = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/processed/longranger_out"
	output:
		filename = config["parent_dir"] + "/data/" + config["exp"] + "_" + config["platform"] + "/processed/longranger_clean/" + config["platform"] + ".R1.gems.no_umi.fa"
	shell:
		".pipeline/barcodes/B_clean_longranger_output.sh {input.directory} {output.filename} {config[platform]}"



