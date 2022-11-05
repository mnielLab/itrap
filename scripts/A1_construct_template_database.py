#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

transtable = str.maketrans("ACTG","TGAC")

"""
Construct a 1OS database of all possible reads from the different combinations of oligo A, oligo B, primers, and samples
"""


# # Input files
OLIGO_1OS = snakemake.input.oligo_1os
OLIGO_CDX = snakemake.input.oligo_cdx
OLIGO_MHC = snakemake.input.oligo_mhc

#OVERREP_SEQS = "../data/" + EXP + "/fastqc/IONTORRENT.R1.gems.no_umi.no_adapters_fastqc/overrepresented_sequences.txt"
NOMENCLATURE = snakemake.input.nomencl

# ## Output files
BARCODE_TEMPLATES = snakemake.output[0] #"../data/" + EXP + "/barcode_library/barcode_templates.fa"

# # Import data
os1 = SeqIO.to_dict(SeqIO.parse(OLIGO_1OS, "fasta"))
mhc = SeqIO.to_dict(SeqIO.parse(OLIGO_MHC, "fasta"))
cdx = SeqIO.to_dict(SeqIO.parse(OLIGO_CDX, "fasta"))

mhc_nom = pd.read_excel(NOMENCLATURE, usecols=[0], names=['Barcode'], sheet_name=0)
cdx_nom = pd.read_excel(NOMENCLATURE, usecols=[0], names=['Barcode'], sheet_name=1)
relevant_barcodes = set(mhc_nom.Barcode.to_list() + cdx_nom.Barcode.dropna().to_list()) # Add CDX barcodes to here


print(mhc)
# # Construct templates
# We will exclude sample IDs since they pollute the annotations
template_records = list()
for barcode_id in mhc_nom.Barcode:
    template_seq = (Seq('N' * 16) +
        Seq('N' * 10) +
        os1['Capture_seq'].seq +
        Seq('N' * 9) +
        mhc[barcode_id].seq +
        Seq('N' * 10) +
        os1['Partial_read_2N'].seq +
        Seq('N' * 10) +
        os1['P7'].seq)
    template_records.append(SeqRecord(template_seq, id=barcode_id))
    print(SeqRecord(template_seq, id=barcode_id).format("fasta"))

for barcode_id in cdx_nom.Barcode:
    template_seq = (Seq('N' * 16) +
        Seq('N' * 10) +
        os1['Capture_seq'].seq +
        Seq('N' * 9) +
        cdx[barcode_id].seq +
        Seq('N' * 10) +
        os1['Partial_read_2N'].seq +
        Seq('N' * 10) +
        os1['P7'].seq)
    template_records.append(SeqRecord(template_seq, id=barcode_id))
    print(SeqRecord(template_seq, id=barcode_id).format("fasta"))


SeqIO.write(template_records, BARCODE_TEMPLATES, "fasta") #Reads are in forward and reverse sense. Where do I find the barcodes? R2?
#SeqIO.write(template_records, EXPECTED_TEMPLATES, "fasta")


