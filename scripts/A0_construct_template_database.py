#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

transtable = str.maketrans("ACTG","TGAC")

"""
Construct a 2OS database of all possible reads from the different combinations of oligo A, oligo B, primers, and samples
"""


# # Args
# TSO, PRIMER_R4, B_OLIGO are same sense as iontorrent seq data
# ANNEAL, A_OLIGO, PRIMER_P1 are negative sense of seq data
N6 = "NNNNNN"
TSO = "TTTCTTATATGGG"
PRIMER_R4 = "CGAGTACCATGGGCGTAAC"
ANNEAL = "GTGTGACCTTCCCCAAAAGGCGTAG".translate(transtable)[::-1]
PRIMER_P1_MHC = "GAAGTTCCAGCCAGCGTC".translate(transtable)[::-1] # Used only for MHC barcodes
PRIMER_P1_CD8 = "GTAAAAGATCCCAGGTTTCATC".translate(transtable)[::-1] # Used only for CD8 barcodes
CD8_BARCODE = "A4000"


# # Input files
OLIGO_A_SEQS = snakemake.input.oligo_a #"../data/" + EXP + "/barcode_library/oligo_a.fa"
OLIGO_B_SEQS = snakemake.input.oligo_b #"../data/" + EXP + "/barcode_library/oligo_b.fa"
SAMPLE_SEQS = snakemake.input.samples #"../data/" + EXP + "/barcode_library/sample.fa"
#OVERREP_SEQS = "../data/" + EXP + "/fastqc/IONTORRENT.R1.gems.no_umi.no_adapters_fastqc/overrepresented_sequences.txt"
NOMENCLATURE = snakemake.input.nomencl

# ## Output files
BARCODE_TEMPLATES = snakemake.output[0] #"../data/" + EXP + "/barcode_library/barcode_templates.fa"

#EXPECTED_TEMPLATES = "../data/" + EXP + "/blast/expected_templates/templates.fa"
#REVERSED_TEMPLATES = "../data/" + EXP + "/blast/reversed_templates/templates.fa"
#REV_TEMPLATES_START = "../data/" + EXP + "/blast/rev_templates_start/templates.fa"
#OVERREP_SEQ_TEMPLATES = "../data/" + EXP + "/blast/overrep_seq_templates/templates.fa"

# # Import data
oligo_a_records = list()
for record in SeqIO.parse(OLIGO_A_SEQS, "fasta"):
    record.seq = record.seq.reverse_complement()
    oligo_a_records.append(record)

oligo_b_records = list()
for record in SeqIO.parse(OLIGO_B_SEQS, "fasta"):
    oligo_b_records.append(record)

sample_records = list()
for record in SeqIO.parse(SAMPLE_SEQS, "fasta"):
    record.seq = record.seq.reverse_complement()
    sample_records.append(record)

mhc_nom = pd.read_excel(NOMENCLATURE, usecols=[0], names=['Barcode'], sheet_name=0)
cdx_nom = pd.read_excel(NOMENCLATURE, usecols=[0], names=['Barcode'], sheet_name=1)
relevant_barcodes = set(mhc_nom.Barcode.to_list() + cdx_nom.Barcode.dropna().to_list()) # Add CDX barcodes to here
#assert 'A4000B288' in relevant_barcodes


# # Construct templates
# We will exclude sample IDs since they pollute the annotations
template_records = list()
for oligo_a_record in oligo_a_records:
    if oligo_a_record.id.endswith(CD8_BARCODE):
        primer_p1 = PRIMER_P1_CD8
    else:
        primer_p1 = PRIMER_P1_MHC           
    for oligo_b_record in oligo_b_records:
        #for sample_record in sample_records:
        template_id = oligo_a_record.id.split("-")[-1] + oligo_b_record.id.split("-")[-1] #+ "_" + sample_record.id
        if template_id in relevant_barcodes:
            template_seq = Seq(TSO + PRIMER_R4 + N6) + oligo_b_record.seq + Seq(ANNEAL) + oligo_a_record.seq + Seq(N6 + primer_p1) #+ sample_record.seq
            template_records.append(SeqRecord(template_seq, id=template_id))
            print(SeqRecord(template_seq, id=template_id).format("fasta"))


SeqIO.write(template_records, BARCODE_TEMPLATES, "fasta")
#SeqIO.write(template_records, EXPECTED_TEMPLATES, "fasta")


