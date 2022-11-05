import glob
import os
import re
import pandas as pd
from scipy.io import mmwrite

configfile: "config/config.13.yaml"
config['run'] = 'run1'
print(config['run'])

include: '/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run1/lib/fq_files.py'
MHC_DIRECTORY = '/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run1/mhc' # change to brc
LIB_DIRECTORY = '/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp13/run1/lib'

rule all:
	input:
		expand(MHC_DIRECTORY + "/mapping/{custom_barcode}_{sorting}/mapping.frag.gz", sorting=sorting.keys(), custom_barcode=custom_barcodes),
		#MHC_DIRECTORY + "/mapping/umi.tsv",
        #MHC_DIRECTORY + "/mapping_split/align/aggregated.aln",
        MHC_DIRECTORY + "/count/matrix.mtx"

def get_filenames(wildcards):
	files = dict()
	for read in ['R1','R2']:
		files[read] = list()
		for f, relevant in zip(['mhc','hsh','mrk'], [mhc_custom, hsh_custom, mrk_custom]):
			if relevant:
				name = sorting[wildcards.sorting][f]

				files[read] += glob.glob(f'{dir_fqs}/**/{name}*{read}*.gz')
	return files

#rule merge_fqs:
#	input:
#		unpack(get_filenames)
#	output:
#		R1 = MHC_DIRECTORY + '/brc_reads.{sorting}.R1.fq.gz',
#		R2 = MHC_DIRECTORY + '/brc_reads.{sorting}.R2.fq.gz'
#	run:
#		#if not input.R2:
#		#	# make R2 if it doesnt exist
#		## merge files
#		shell("ls {input}")
#		if len(input.R1) == len(input.R2):
#			shell("cat {input.R1} > {output.R1}")
#			shell("cat {input.R2} > {output.R2}")
#		if len(input.R1) != len(input.R2):
#			regex = '_S\d+_L\d+_R\d+_\d+\.\w+\.gz'
#
#			R1_files = set()
#			for filename in input.R1:
#				R1_files.add(re.sub(regex, '', os.path.basename(filename)))
#
#			R2_files = set()
#			for filename in input.R2:
#				R2_files.add(re.sub(regex, '', os.path.basename(filename)))
#
#			R2_missing = R1_files - R2_files
#
#			for filename in input.R1:
#				if any(miss in filename for miss in R2_missing):
#					# make dummy read

#################################################################
#                   Barcodes analyzed by KMA                    #
#################################################################

def get_filenames(wildcards):
	files = dict()
	for read in ['R1','R2']:
		prefix = sorting[wildcards.sorting][wildcards.custom_barcode]
		file_lst = glob.glob(f'{dir_fqs}/**/{prefix}*{read}*.gz')
		assert len(file_lst) == 1
		files[read] = file_lst[0]
	return files

rule cut_adapters:
	"""
	If the sequencing provider is IONTORRENT for the custom barcodes
	it is important to trim the 5' adapters before running longranger
	to fetch the GEM barcodes.
	I assume that if the data is single end then it is also IONTORRENT.
	If the data is single end I also need to write a dummy R2 for longranger.
	"""
	input:
		unpack(get_filenames)
	output:
		R1 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
		R2 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2'])
	params:
		ADAPTER_5p = "CACGACGCTCTTCCGATCT", # IONTORRENT
		cutadapt = config['cutadapt']
	run:
		print(input)
		if input.R2:
			# Paired end reads
			shell("ln -sr {input.R1} {output.R1}")
			shell("ln -sr {input.R2} {output.R2}")
		else:
			# Single end reads
			shell("{params.cutadapt} -g {params.ADAPTER_5p} -o {output.R1} {input.R1}")
			# Make dummy read
			shell("./workflow/wrappers/wrt_dummy_reads.sh {input.R1} {output.R2}")
			# zcat $1 | awk --posix '{ if (NR % 4 == 0) { sub(/.*/,"CCCCCCCC")} else if (NR % 2 == 0) { sub(/.*/,"NNNNNNNN")}; print}' > $2

rule run_longranger:
	# Does not handle single end reads, only paired end.
	input:
		R1 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
		R2 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2'])
	output:
		touch(MHC_DIRECTORY + '/longranger/{custom_barcode}_{sorting}.done')
	params:
		rundir = lambda wildcards, output: os.path.dirname(os.path.abspath(output[0])),
		longranger = config['longranger'],
		identifier = '{custom_barcode}_{sorting}',
		fastqs = lambda wildcards, input: os.path.dirname(os.path.abspath(input.R1)), #dir_fqs,
		sample = lambda wildcards, input: re.sub('_S\d+_L\d+_R\d+_\d+\.\w+\.gz', '', os.path.basename(input.R1)) #sorting[wildcards.sorting][wildcards.custom_barcode]
	shell:
		"""
		cd {params.rundir}

		{params.longranger} basic --id {params.identifier} --fastqs {params.fastqs} --sample {params.sample}
		"""

rule deinterleave_longranger_output:
	# Should I also remove the remaining 3 UMI?
	input:
		MHC_DIRECTORY + '/longranger/{custom_barcode}_{sorting}.done'
	params:
		interleaved_fq = MHC_DIRECTORY + '/longranger/{custom_barcode}_{sorting}/outs/barcoded.fastq.gz'
	output:
		R1 = MHC_DIRECTORY + '/longranger/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
		R2 = MHC_DIRECTORY + '/longranger/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2'])
	shell:
		"""
		./workflow/wrappers/deinterleave_fastq.sh {params.interleaved_fq} {output}
		"""

rule make_templates:
	input:
		brc_info = os.path.join(LIB_DIRECTORY, 'barcode_design.yaml'),
		nomencl = os.path.join(LIB_DIRECTORY, 'barcode_specificity_annotations.xlsx')
	output:
		MHC_DIRECTORY + "/mapping/barcode_templates.fa"
	conda:
		"../envs/biopython.yaml"
	script:
		"../../scripts/A2_construct_template_database.py"

rule run_kma_index:
	input:
		MHC_DIRECTORY + "/mapping/barcode_templates.fa"
	output:
		MHC_DIRECTORY + "/mapping/templates.name"
	params:
		prefix = lambda wildcards, output: os.path.splitext(output[0])[0],
		mapper = config['KMA']
	shell:
		"{params.mapper} index -i {input} -o {params.prefix}"

rule run_kma_map:
	input:
		R1 = MHC_DIRECTORY + '/longranger/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
		R2 = MHC_DIRECTORY + '/longranger/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2']),
		templates = MHC_DIRECTORY + "/mapping/templates.name"
	output:
		MHC_DIRECTORY + "/mapping/{custom_barcode}_{sorting}/mapping.frag.gz"
	params:
		input_prefix = lambda wildcards, input: os.path.splitext(input.templates)[0], #[-1]
		output_prefix = lambda wildcards, output: output[0].rsplit('.', 2)[0],
		mapper = config['KMA'],
		arguments = config["mapping_args"]
	shell:
		"""
		set +o pipefail;

		dummy_R2=$(zcat {input.R2} | grep 'NNNNNNNN' | wc -l) #| paste - - - - 

		if [ $dummy_R2 -eq 0 ]; then
			echo Paired end reads
			{params.mapper} -ipe {input.R1} {input.R2} -o {params.output_prefix} -t_db {params.input_prefix} {params.arguments}
		else
			echo Single end reads
			{params.mapper} -i {input.R1} -o {params.output_prefix} -t_db {params.input_prefix} {params.arguments}
		fi
		"""
rule clean_and_merge_kma_output:
    # Maybe still parallelize on custom_barcode?
	input:
		expand(MHC_DIRECTORY + "/mapping/{custom_barcode}_{sorting}/mapping.frag.gz", sorting=sorting.keys(), custom_barcode=custom_barcodes)
	output:
		MHC_DIRECTORY + "/mapping/mapping.frag.gz"
	run:
		dfs = list()
		for filename in input:
			matches = re.search('/mapping/(\w+)_(\w+)/mapping.frag.gz', filename)
			barcode = matches.group(1)
			sorting = matches.group(2)

			try:
				tmp = pd.read_csv(filename, sep='\t', header=None)
				tmp['brc'] = barcode
				tmp['label'] = 1 if 'pos' in sorting else 1 if 'tot' in sorting else 0

				# select only unambiguous hits
				dfs.append(tmp[tmp[1] <= 1])

			except pd.errors.EmptyDataError:
				print(f'{filename} is empty')

		print(output[0])
		pd.concat(dfs, ignore_index=True).to_csv(output[0], sep='\t', header=False, index=False)

# --cores 20
rule extract_umis:
	input:
		MHC_DIRECTORY + "/mapping/mapping.frag.gz",
		"tools/barcode_oligos/samples.rv.fa", #config['samples']['default'], #MHC_DIRECTORY + "/barcode_library/sample.fa",
		"tools/barcode_oligos/oligo_a.25mer.rv.fa", #config['oligo_a']['25mer'], #MHC_DIRECTORY + "/barcode_library/oligo_a.fa",
		config['oligo_b']['25mer'], #MHC_DIRECTORY + "/barcode_library/oligo_b.fa",
		"tools/barcode_oligos/barcode-information.fa"
	output:
		MHC_DIRECTORY + "/mapping/umi.tsv"
	#conda:
	#	"snake_envs/umi.yaml"
	shell:
		"/home/tuba-nobackup/shared/R/R-3.6.1/bin/Rscript ./scripts/parse-kma-results.R {input} {output}"
#

#rule detect_umis:
#
#rule augment_barcode_data:

# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint split_mapping_output_by_template:
    input:
        MHC_DIRECTORY + "/mapping/mapping.frag.gz"
    output:
        directory(MHC_DIRECTORY + "/mapping/split/map/")
    shell:
        """
        zcat {input} | awk '{{print>"{output}"$6}}'
        """

# an intermediate rule
rule prep_pairwise_alignment_inputs:
    input:
        aseq = MHC_DIRECTORY + "/mapping/templates.fa",
        bseq = MHC_DIRECTORY + "/mapping/split/map/{i}" #map
    output:
        aseq = MHC_DIRECTORY + "/mapping/split/seq/{i}.aseq.fa", #seq
        bseq = MHC_DIRECTORY + "/mapping/split/seq/{i}.bseq.fa"
    params:
        i = '{i}'
    shell:
        """
        grep '{params.i}' -A 1 {input.aseq} > {output.aseq}
 		awk -F'\t' -v OFS='\t' '{{print $7, $1}}' {input.bseq} | awk -F'\t' -v OFS='\n' '{{$1 = ">" $1}} 1' > {output.bseq}
 		"""
        
rule pairwise_alignment:
    input:
        aseq = MHC_DIRECTORY + "/mapping/split/seq/{i}.aseq.fa",
        bseq = MHC_DIRECTORY + "/mapping/split/seq/{i}.bseq.fa"
    output:
        outf = MHC_DIRECTORY + "/mapping/split/aln/{i}.aln" #aln
    params:
        aligner = config['pairwise_aligner'],
        o = 10,
        e = 0.5
    shell:
        """
        {params.aligner} -asequence {input.aseq} -bsequence {input.bseq} -outfile {output.outf} -gapopen {params.o} -gapextend {params.e}
        """

rule extract_umi:
    input:
        MHC_DIRECTORY + "/mapping/split/aln/{i}.aln"
    output:
        MHC_DIRECTORY + "/mapping/split/umi/{i}.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../../scripts/extract_umi.py"
        
rule count_umi:
    input:
        mpp = MHC_DIRECTORY + "/mapping/split/map/{i}",
        umi = MHC_DIRECTORY + "/mapping/split/umi/{i}.tsv"
    output:
        MHC_DIRECTORY + "/mapping/split/cnt/{i}.csv"
    script:
        "../../scripts/count_umi.py"
        
#def aggregate_input(wildcards):
#    checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
#    return expand(MHC_DIRECTORY + "/mapping_split/{d}/{i}",
#           d=['templ','mappd'],
#           i=glob_wildcards(os.path.join(checkpoint_output, "{i}")).i)

#def aggregate_input(wildcards):
#    checkpoint_output = checkpoints.split_mapping_output_by_template.get(**wildcards).output[0]
#    meh = expand(MHC_DIRECTORY + "/mapping_split/mappd/{i}.bseq.fa", i=glob_wildcards(os.path.join(checkpoint_output, "{i}")).i) #"{i, (?!snakemake_timestamp).*}.bseq.fa"
#    print(meh)
#    print([m for m in meh if re.match('^((?!snakemake_timestamp).)*$', m)])
#    return [m for m in meh if re.match('^((?!snakemake_timestamp).)*$', m)]

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_mapping_output_by_template.get(**wildcards).output[0]
    template_ids = glob_wildcards(os.path.join(checkpoint_output, "{i}")).i
    # For some reason the function picks up the timestamp files?!
    i = [t for t in template_ids if re.match('^((?!snakemake_timestamp).)*$', t)]
    print(i)
    return expand(MHC_DIRECTORY + "/mapping/split/cnt/{i}.csv", i=i)

# an aggregation over all produced clusters
rule aggregate:
    input:
        aggregate_input
    output:
        temp(MHC_DIRECTORY + "/mapping/split/umi.tsv")
    shell:
        "cat {input} > {output.tmp}"
        
rule write_count_matrix:
    input:
        MHC_DIRECTORY + "/mapping/split/umi.tsv"
    output:
        mtx = MHC_DIRECTORY + "/count/matrix.mtx",
        brc = MHC_DIRECTORY + "/count/features.tsv",
        gem = MHC_DIRECTORY + "/count/barcodes.tsv"
    script:
        "../../scripts/wrt_count_matrix.py"
        
rule parse_count_matrix:
    # Take input from either cellranger or KMA.
    input:
        mtx = MHC_DIRECTORY + "/count/matrix.mtx",
        brc = MHC_DIRECTORY + "/count/features.tsv",
        gem = MHC_DIRECTORY + "/count/barcodes.tsv"
    params:
        mpp = MHC_DIRECTORY + "/mapping/mapping.frag.gz"
    output:
        
    script:
        "../../scripts/prs_count_matrix.py"
        
# How to compile output from KMA and cellranger?
# Use rule inheritance:
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-inheritance
# Use Handling Ambiguous Rules
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#handling-ambiguous-rules


