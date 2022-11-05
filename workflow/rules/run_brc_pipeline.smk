import pandas as pd
import re
import os
#################################################################
#                            Targets                            #
#################################################################
TARGET['write_count_matrix'] = MHC_DIRECTORY + "/count/matrix.kma.mtx"

#################################################################
#                             Rules                             #
#################################################################

#################################################################
#                   Barcodes analyzed by KMA                    #
#################################################################

#def get_filenames(wildcards):
#    files = dict()
#    for read in ['R1','R2']:
#        prefix = sorting[wildcards.sorting][wildcards.custom_barcode]
#        print(prefix)
#        print(f'{dir_fqs}/**/{prefix}*{read}*.gz')
#        file_lst = glob.glob(f'{dir_fqs}/**/{prefix}*{read}*.gz')
#        print(file_lst)
#        if file_lst == []:
#            print(glob.glob(f'{dir_fqs}/{prefix}*.gz'))
#            file_lst = glob.glob(f'{dir_fqs}/{prefix}*{read}*.gz')
#        print(file_lst)
#        assert len(file_lst) == 1
#        files[read] = file_lst[0]
#    return files

def get_filenames(wildcards):
    prefix = sorting[wildcards.sorting][wildcards.custom_barcode]
    file_lst = glob.glob(f'{dir_fqs}/**/{prefix}*.gz', recursive=True)
    
    files = dict()
    for filename in file_lst:
        read = re.search('_S\d+_L\d+_(R\d+)_\d+\.\w+\.gz', filename).group(1)
        files[read] = filename
        
    print(files)
    return files
    
ruleorder: cut_adapters > skip_adapter_trimming
#rule cut_adapters:
#    """
#    If the sequencing provider is IONTORRENT for the custom barcodes
#    it is important to trim the 5' adapters before running longranger
#    to fetch the GEM barcodes.
#    I assume that if the data is single end then it is also IONTORRENT.
#    If the data is single end I also need to write a dummy R2 for longranger.
#    """
#    input:
#        unpack(get_filenames)
#    output:
#        R1 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
#        R2 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2'])
#    params:
#        ADAPTER_5p = "CACGACGCTCTTCCGATCT", # IONTORRENT
#        cutadapt = config['cutadapt']
#    run:
#        print('running cut_adapters')
#        if input.R2:
#            print('linking the paired reads')
#            # Paired end reads
#            shell("ln -sr {input.R1} {output.R1}")
#            shell("ln -sr {input.R2} {output.R2}")
#        else:
#            print('making dummy R2 reads')
#            # Single end reads
#            shell("{params.cutadapt} -g {params.ADAPTER_5p} -o {output.R1} {input.R1}")
#            # Make dummy read
#            shell("./workflow/wrappers/wrt_dummy_reads.sh {input.R1} {output.R2}")
#            # zcat $1 | awk --posix '{ if (NR % 4 == 0) { sub(/.*/,"CCCCCCCC")} else if (NR % 2 == 0) { sub(/.*/,"NNNNNNNN")}; print}' > $2

rule cut_adapters:
    input:
        unpack(get_filenames)
    output:
        R1 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
        R2 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2'])
    params:
        ADAPTER_5p = "CACGACGCTCTTCCGATCT", # IONTORRENT
        cutadapt = config['cutadapt']
    run:
        print('making dummy R2 reads')
        # Single end reads
        shell("{params.cutadapt} -g {params.ADAPTER_5p} -o {output.R1} {input.R1}")
        # Make dummy read
        shell("./workflow/wrappers/wrt_dummy_reads.sh {input.R1} {output.R2}")
        # zcat $1 | awk --posix '{ if (NR % 4 == 0) { sub(/.*/,"CCCCCCCC")} else if (NR % 2 == 0) { sub(/.*/,"NNNNNNNN")}; print}' > $2
        
rule skip_adapter_trimming:
    input:
        unpack(get_filenames)
    output:
        R1 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
        R2 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2'])
    run:
        print('Skipping adapter trimming and linking the paired reads')
        # Paired end reads
        shell("ln -sr {input.R1} {output.R1}")
        shell("ln -sr {input.R2} {output.R2}")
            
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
    # Maybe just merge without annotating custom barcode and sorting?
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

## --cores 20
#rule extract_umis:
#	input:
#		MHC_DIRECTORY + "/mapping/mapping.frag.gz",
#		"tools/barcode_oligos/samples.rv.fa", #config['samples']['default'], #MHC_DIRECTORY + "/barcode_library/sample.fa",
#		"tools/barcode_oligos/oligo_a.25mer.rv.fa", #config['oligo_a']['25mer'], #MHC_DIRECTORY + "/barcode_library/oligo_a.fa",
#		config['oligo_b']['25mer'], #MHC_DIRECTORY + "/barcode_library/oligo_b.fa",
#		"tools/barcode_oligos/barcode-information.fa"
#	output:
#		MHC_DIRECTORY + "/mapping/umi.tsv"
#	#conda:
#	#	"snake_envs/umi.yaml"
#	shell:
#		"/home/tuba-nobackup/shared/R/R-3.6.1/bin/Rscript ./scripts/parse-kma-results.R {input} {output}"
#

#rule detect_umis:
#
#rule augment_barcode_data:

# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint split_mapping_output_by_template:
    input:
        MHC_DIRECTORY + "/mapping/mapping.frag.gz"
    output:
        directory(MHC_DIRECTORY + "/mapping/spl/") #split/map/
    shell:
        """
        # Many errors! This step might have to be done manually :(
        #set +e
        #set +o pipefail
        cd {output}
        echo $?
        echo $(ls {output})
        echo $?
        zcat {input} | awk '{{print > $6}}'
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """
        #./workflow/wrappers/split_fasta.sh {input} {output}
#        if [ ! -d {output} ]; then mkdir {output}; fi

# an intermediate rule
rule prep_pairwise_alignment_inputs:
    input:
        aseq = MHC_DIRECTORY + "/mapping/barcode_templates.fa",
        bseq = MHC_DIRECTORY + "/mapping/spl/{i}" #map
    output:
        aseq = MHC_DIRECTORY + "/mapping/seq/{i}.aseq.fa", #seq
        bseq = MHC_DIRECTORY + "/mapping/seq/{i}.bseq.fa"
    params:
        i = '{i}'
    shell:
        """
        grep '{params.i}' -A 1 {input.aseq} > {output.aseq}
        awk -F'\t' -v OFS='\t' '{{print $7, $1}}' {input.bseq} | awk -F'\t' -v OFS='\n' '{{$1 = ">" $1}} 1' > {output.bseq}
        """
        
rule pairwise_alignment:
    input:
        aseq = MHC_DIRECTORY + "/mapping/seq/{i}.aseq.fa", #seq
        bseq = MHC_DIRECTORY + "/mapping/seq/{i}.bseq.fa"
    output:
        outf = MHC_DIRECTORY + "/mapping/aln/{i}.aln" #aln
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
        MHC_DIRECTORY + "/mapping/aln/{i}.aln"
    output:
        MHC_DIRECTORY + "/mapping/umi/{i}.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../../scripts/extract_umi.py"
        
rule count_umi:
    input:
        mpp = MHC_DIRECTORY + "/mapping/spl/{i}",
        umi = MHC_DIRECTORY + "/mapping/umi/{i}.tsv"
    output:
        MHC_DIRECTORY + "/mapping/cnt/{i}.csv"
    script:
        "../../scripts/count_umi.py"
        
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_mapping_output_by_template.get(**wildcards).output[0]
    print(checkpoint_output)
    template_ids = glob_wildcards(os.path.join(checkpoint_output, "{i}")).i
    # For some reason the function picks up the timestamp files?!
    i = [t for t in template_ids if re.match('^((?!snakemake_timestamp).)*$', t)]
    print('aggregate func')
    print(i)
    #print(expand(MHC_DIRECTORY + "/mapping/cnt/{i}.csv", i=i))
    return expand(MHC_DIRECTORY + "/mapping/cnt/{i}.csv", i=i)

tmp_barcode_files = ['A4000B303', 'A1072B303', 'A1073B303', 'A1074B303', 'A1075B303', 'A1076B303', 'A1077B303', 'A1078B303', 'A1079B303', 'A1072B304', 'A1073B304', 'A1074B304', 'A1075B304', 'A1076B304', 'A1077B304', 'A1078B304', 'A1072B305', 'A1073B305', 'A1074B305', 'A1075B305', 'A1076B305', 'A1077B305', 'A1078B305', 'A1079B305', 'A1072B306', 'A1073B306', 'A1074B306', 'A1075B306', 'A1076B306', 'A1077B306', 'A1078B306', 'A1079B306', 'A1072B307', 'A1073B307', 'A1074B307', 'A1075B307', 'A1076B307', 'A1077B307', 'A1078B307', 'A1079B307', 'A1072B308', 'A1073B308', 'A1075B308', 'A1077B308', 'A1078B308', 'A1079B308', 'A4000B309', 'A1072B309', 'A1073B309', 'A1074B309', 'A1075B309', 'A1076B309', 'A1077B309', 'A1078B309', 'A1079B309', 'A1072B310', 'A1073B310', 'A1074B310', 'A1075B310', 'A1076B310', 'A1077B310', 'A1078B310', 'A1079B310', 'A4000B311', 'A1072B311', 'A1073B311', 'A1074B311', 'A1075B311', 'A1076B311', 'A1077B311', 'A1078B311', 'A1079B311', 'A1072B313', 'A1073B313', 'A1074B313', 'A1075B313', 'A1076B313', 'A1077B313', 'A1078B313', 'A1079B313', 'A1072B314', 'A1073B314', 'A1074B314', 'A1075B314', 'A1076B314', 'A1077B314', 'A1078B314', 'A1079B314', 'A4000B315', 'A1072B315', 'A1073B315', 'A1074B315', 'A1075B315', 'A1076B315', 'A1077B315', 'A1078B315', 'A1079B315', 'A4000B316', 'A1072B316', 'A1073B316', 'A1074B316', 'A1076B316', 'A1077B316', 'A1079B316', 'A1072B317', 'A1073B317', 'A1074B317', 'A1075B317', 'A1077B317', 'A1078B317', 'A1072B318', 'A1073B318', 'A1074B318', 'A1075B318', 'A1076B318', 'A1077B318', 'A1078B318', 'A1079B318', 'A4000B304', 'A4000B305', 'A4000B306', 'A4000B307', 'A4000B310', 'A4000B313', 'A4000B314', 'A4000B317', 'A4000B318']

# an aggregation over all produced clusters
rule aggregate_alignments:
    input:
        #aggregate_input # <- problems with this step! Manually get a list of the barcodes
        expand(MHC_DIRECTORY + "/mapping/cnt/{i}.csv", i=tmp_barcode_files)
    output:
        MHC_DIRECTORY + "/mapping/umi.tsv" # temp()
    shell:
        "cat {input} > {output}"
        
rule write_count_matrix:
    input:
        MHC_DIRECTORY + "/mapping/umi.tsv"
    output:
        brc = MHC_DIRECTORY + "/count/features.kma.tsv",
        gem = MHC_DIRECTORY + "/count/barcodes.kma.tsv",
        mtx = MHC_DIRECTORY + "/count/matrix.kma.mtx"
    script:
        "../../scripts/wrt_count_matrix.py"
        
