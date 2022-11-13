import pandas as pd
import re
import os

#################################################################
#                            Targets                            #
#################################################################
"""
Add intermediate targets to pipeline to test partial run, e.g. add:

TARGET['write_count_matrix'] = MHC_DIRECTORY + "/count/matrix.kma.mtx"
"""

TARGET['write_count_matrix'] = MHC_DIRECTORY + "/count/matrix.kma.mtx"

#################################################################
#                             Rules                             #
#################################################################

#################################################################
#                   Barcodes analyzed by KMA                    #
#################################################################

def get_filenames(wildcards):
    """
    Retrieving files based on prefixes listed in config.
    """
    prefix = sorting[wildcards.sorting][wildcards.custom_barcode]
    file_lst = glob.glob(f'{dir_fqs}/**/{prefix}*.gz', recursive=True)
    
    files = dict()
    for filename in file_lst:
        read = re.search('_S\d+_L\d+_(R\d+)_\d+\.\w+\.gz', filename).group(1)
        files[read] = filename
        
    return files
    
ruleorder: cut_adapters > skip_adapter_trimming

rule cut_adapters:
    """
    Trimming adapters using cutadapt.
    Assuming barcode reads are not paired end.
    Making dummy reads for R2, since Longranger does not handle single-end reads, but only paired end.
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
        print('making dummy R2 reads')
        # Single end reads
        shell("{params.cutadapt} -g {params.ADAPTER_5p} -o {output.R1} {input.R1}")
        # Make dummy read
        shell("./workflow/wrappers/wrt_dummy_reads.sh {input.R1} {output.R2}")
        
rule skip_adapter_trimming:
    """
    Skipping adapter trimming.
    Assuming paired end reads.
    """
    input:
        unpack(get_filenames)
    output:
        R1 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
        R2 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2'])
    run:
        shell("ln -sr {input.R1} {output.R1}")
        shell("ln -sr {input.R2} {output.R2}")
            
rule run_longranger:
    """
    Running Longranger on R1 and (dummy) R2 to extract 10x barcode.
    Longranger partially removes the UMI directly upstream of the 10x barcode.
    Longranger produces output dir of random name in CWD.
    """
    input:
        R1 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R1']),
        R2 = MHC_DIRECTORY + '/cut_adapters/{{custom_barcode}}_{{sorting}}_{ext}'.format(ext=config['dummy_ext']['R2'])
    output:
        touch(MHC_DIRECTORY + '/longranger/{custom_barcode}_{sorting}.done')
    params:
        rundir = lambda wildcards, output: os.path.dirname(os.path.abspath(output[0])),
        longranger = config['longranger'],
        identifier = '{custom_barcode}_{sorting}',
        fastqs = lambda wildcards, input: os.path.dirname(os.path.abspath(input.R1)),
        sample = lambda wildcards, input: re.sub('_S\d+_L\d+_R\d+_\d+\.\w+\.gz', '', os.path.basename(input.R1))
    shell:
        """
        cd {params.rundir}
        
        {params.longranger} basic --id {params.identifier} --fastqs {params.fastqs} --sample {params.sample}
        """

rule deinterleave_longranger_output:
    """
    Longranger outputs interleaved R1 and R2 reads which are separated into two files again.
    The remaining UMI are removed in this step.
    """
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
    """
    Construct barcode reference sequences comprising of the DNA barcode library.
    """
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
    """
    Index the barcode reference sequences for fast mapping.
    """
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
    """
    Map DNA barcode reads against barcode reference sequences.
    """
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
    """
    Merge barcodes, when both positive and negative staining pools are sequenced.
    Remove reads mapping to multiple reference templates.
    """
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


# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint split_mapping_output_by_template:
    """
    Split file of mapped barcodes into file per barcode.
    """
    input:
        MHC_DIRECTORY + "/mapping/mapping.frag.gz"
    output:
        directory(MHC_DIRECTORY + "/mapping/spl/")
    shell:
        """
        # Many errors! This step might have to be done manually :(
        #set +e
        #set +o pipefail
        cd {output}
        zcat {input} | awk '{{print > $6}}'
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """

# an intermediate rule
rule prep_pairwise_alignment_inputs:
    """
    Split barcode reference fasta into a fasta file per template.
    Convert table of mapped results to fasta format per barcode.
    """
    input:
        aseq = MHC_DIRECTORY + "/mapping/barcode_templates.fa",
        bseq = MHC_DIRECTORY + "/mapping/spl/{bc}"
    output:
        aseq = MHC_DIRECTORY + "/mapping/seq/{bc}.aseq.fa",
        bseq = MHC_DIRECTORY + "/mapping/seq/{bc}.bseq.fa"
    params:
        bc = '{bc}'
    shell:
        """
        grep '{params.bc}' -A 1 {input.aseq} > {output.aseq}
        awk -F'\t' -v OFS='\t' '{{print $7, $1}}' {input.bseq} | awk -F'\t' -v OFS='\n' '{{$1 = ">" $1}} 1' > {output.bseq}
        """
        
rule pairwise_alignment:
    """
    Align reads against template to identify UMI sequences and quantify abundance.
    """
    input:
        aseq = MHC_DIRECTORY + "/mapping/seq/{bc}.aseq.fa",
        bseq = MHC_DIRECTORY + "/mapping/seq/{bc}.bseq.fa"
    output:
        outf = MHC_DIRECTORY + "/mapping/aln/{bc}.aln"
    params:
        aligner = config['pairwise_aligner'],
        o = 10,
        e = 0.5
    shell:
        """
        {params.aligner} -asequence {input.aseq} -bsequence {input.bseq} -outfile {output.outf} -gapopen {params.o} -gapextend {params.e}
        """

rule extract_umi:
    """
    Identify and extract UMI to quantify abundance.
    """
    input:
        MHC_DIRECTORY + "/mapping/aln/{bc}.aln"
    output:
        MHC_DIRECTORY + "/mapping/umi/{bc}.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../../scripts/extract_umi.py"
        
rule count_umi:
    """
    Quantify abundance of each barcode by counts of UMI.
    """
    input:
        mpp = MHC_DIRECTORY + "/mapping/spl/{bc}",
        umi = MHC_DIRECTORY + "/mapping/umi/{bc}.tsv"
    output:
        MHC_DIRECTORY + "/mapping/cnt/{bc}.csv"
    script:
        "../../scripts/count_umi.py"
        
def aggregate_input(wildcards):
    """
    Produce list of expected output filenames from count_umi rule.
    OBS! Fragile function. Does not behave well upon reruns!
    All files from checkpoint until aggregate function should be removed for proper rerun.
    """
    checkpoint_output = checkpoints.split_mapping_output_by_template.get(**wildcards).output[0]
    template_ids = glob_wildcards(os.path.join(checkpoint_output, "{bc}")).bc
    return expand(MHC_DIRECTORY + "/mapping/cnt/{bc}.csv", bc=template_ids)


# an aggregation over all produced clusters
rule aggregate_alignments:
    """
    Aggregates all UMI counts per barcode per GEM.
    """
    input:
        aggregate_input
        # If problems with above, get a list of barcodes (barcode_file_list) and run line below instead.
        # expand(MHC_DIRECTORY + "/mapping/cnt/{bc}.csv", i=barcode_file_list)
    output:
        MHC_DIRECTORY + "/mapping/umi.tsv"
    shell:
        "cat {input} > {output}"
        
rule write_count_matrix:
    """
    Convert barcode data into a count matrix of row barcodes and columns GEMs.
    """
    input:
        MHC_DIRECTORY + "/mapping/umi.tsv"
    output:
        brc = MHC_DIRECTORY + "/count/features.kma.tsv",
        gem = MHC_DIRECTORY + "/count/barcodes.kma.tsv",
        mtx = MHC_DIRECTORY + "/count/matrix.kma.mtx"
    script:
        "../../scripts/wrt_count_matrix.py"
        
