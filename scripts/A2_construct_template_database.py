#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import yaml
import itertools

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
Construct a any barcode database of all possible reads from the different combinations of oligo A, oligo B, primers, and samples
"""

def load_fasta_headers(filename):
    def get_id_only(identifier):
        return identifier.split('-')[-1]
    return SeqIO.index(filename, 'fasta', key_function=get_id_only)

def get_relevant_oligos(headers):
    all_oligo_ids = pd.concat([mhc_nom.Barcode, cdx_nom.Barcode]).astype(str)
    relevant_oligos = list()
    for header in headers.keys():
        if any(header in oligo for oligo in all_oligo_ids):
            relevant_oligos.append(header)
    return relevant_oligos #set(headers.keys()).intersection(set(all_oligo_ids))

def get_fasta_entries(records, oligo_ids):
    lst = list()
    for oid in oligo_ids:
        lst.append(records[oid])
    return lst

def load_barcode_building_blocks():
    with open(BRC_INFO, 'r') as f:
        return yaml.safe_load(f)

def insert_brc_oligo_dependent_sequences(attributes):
    """
    Some barcode segments may be dependent on the barcode ID, which may reflect if its a pMHC or CDX barcode.
    In one desgin the P1 primer depends on a specific A oligo.
    Due to itertools.product we will generate constructs that are unintended.
    We check for special cases where oligo A and P1 agree.
    Below is listed putative combinations:

    Keep A    B    P1
    T    A4000B759 A4000 dependent
    F    A4234B759 A4000 dependent
    F    A4000B759 default primer
    T    A4234B759 default primer

    There may be multiple barcode segments that contain such a dependency.
    """
    for attribute in attributes: # [P1]
        intended_barcode_design = list()

        special_case = any([tmp_id in bbb[attribute].keys() for tmp_id in template_ids])

        if (template_bbb[attribute] != 'default') & special_case:
            intended_barcode_design.append(True) #Seq(str(template_seq) % oligo_seq) # Insert seq at placeholder
        elif (template_bbb[attribute] == 'default') & (not special_case):
            intended_barcode_design.append(True) #Seq(str(template_seq) % bbb[attribute]['default']) # Insert default seq at placeholder
        else:
            # There is a mismatch between template ID and bbb[attribute].value
            # Not an intended barcode.
            intended_barcode_design.append(False)

    if all(intended_barcode_design):
        return True
    else:
        return False

def product_dict(**kwargs):
    keys = kwargs.keys()
    vals = kwargs.values()
    for instance in itertools.product(*vals):
        yield dict(zip(keys, instance))

# # Input files
#OLIGO_A_SEQS = snakemake.input.oligo_a #"../data/" + EXP + "/barcode_library/oligo_a.fa"
#OLIGO_B_SEQS = snakemake.input.oligo_b #"../data/" + EXP + "/barcode_library/oligo_b.fa"
#SAMPLE_SEQS = snakemake.input.samples #"../data/" + EXP + "/barcode_library/sample.fa"
NOMENCLATURE = snakemake.input.nomencl
BRC_INFO = snakemake.input.brc_info

# ## Output files
BARCODE_TEMPLATES = snakemake.output[0] #"../data/" + EXP + "/barcode_library/barcode_templates.fa"

# Load
mhc_nom = pd.read_excel(NOMENCLATURE, usecols=[0], names=['Barcode'], sheet_name=0)
cdx_nom = pd.read_excel(NOMENCLATURE, usecols=[0], names=['Barcode'], sheet_name=1)

bbb = load_barcode_building_blocks()

# Main

# Prepare barcode building blocks: Compile sequences and convert str to [Seq]
for key, value in bbb.items():
    if key in snakemake.config:
        # Load the relevant contents of the file
        filename = snakemake.config[key][value]
        lib_oligos = load_fasta_headers(filename)
        relevant_oligos = get_relevant_oligos(lib_oligos)
        records = get_fasta_entries(lib_oligos, relevant_oligos)
        # Overwrite filename with actual content
        bbb[key] = records
    elif isinstance(value, str):
        # Convert string to sequence in list
        # Prep for dictionary product expansion
        bbb[key] = [Seq(value)]
    elif isinstance(value, dict):
        # The P1 primer may depend on an oligo and therefore contains a dict
        for k,v in value.items():
            if isinstance(v, str):
                bbb[key][k] = Seq(v)
            else:
                print('ERROR: unexpected format of barcode config')

# Build templates
template_records = list()
for building_blocks in itertools.product(*bbb.values()):
    template_ids = list()
    template_seq = Seq('')
    template_bbb = dict(zip(bbb.keys(), building_blocks))

    attributes, sequences = list(), list()

    for name, bb in template_bbb.items():
        if isinstance(bb, Bio.Seq.Seq):
            template_seq += bb
        elif isinstance(bb, Bio.SeqRecord.SeqRecord):
            template_seq += bb.seq
            template_ids.append(bb.id)
        elif isinstance(bb, str):
            template_seq += bbb[name][bb] # Eg name = P1, bb = A4000
            attributes.append(name) # Eg. P1 can take different values depending on oligo A

            #sequences.append(bbb[name][bb]) # Eg name = P1, bb = A4000
            #template_seq += Seq('%s') # Insert placeholder
            
    print(attributes)

    if insert_brc_oligo_dependent_sequences(attributes):
        #template_seq = Seq(str(template_seq) % tuple(sequences))
        template_records.append(SeqRecord(template_seq, id=''.join(sorted(template_ids))))
    else:
        print('not a match', ''.join(sorted(template_ids)))

#assert len(template_records) >= len(mhc_nom) + len(cdx_nom)
print('N templates:', len(template_records))
print('N barcodes:', len(mhc_nom) + len(cdx_nom))
SeqIO.write(template_records, BARCODE_TEMPLATES, "fasta-2line") # changed to format variant with no line wrapping and exactly two lines per record. # "fasta" 





