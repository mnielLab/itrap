#!/usr/bin/env python
# data/exp9.2_RAW/90-543315613/00_fastq/9thantibody_R1_001.fastq.gz
dir_fqs = "/home/tuba/herpov/tcr-pmhc-sc-project/data/exp9.2_RAW/90-543315613/00_fastq" # OBS! I have removed the trailing /

pos_tcr = "9thTCR"
neg_tcr = ""

pos_mhc = "9thMHCION"
neg_mhc = ""

pos_hsh = "9thHSHION" # Hashing label
neg_hsh = "" # Hashing label

pos_mrk = "" # Surface/feature marker label
neg_mrk = "" # Surface/feature marker label

pos_gex = "9thfulltranscriptome"
neg_gex = ""

mhc_custom = True
hsh_custom = True
mrk_custom = None

if any((neg_tcr, neg_mhc, neg_hsh, neg_mrk, neg_gex)):
    sorting = {'pos': {'tcr': pos_tcr, 'mhc': pos_mhc, 'hsh': pos_hsh, 'mrk': pos_mrk, 'gex': pos_gex},
               'neg': {'tcr': neg_tcr, 'mhc': neg_mhc, 'hsh': neg_hsh, 'mrk': neg_mrk, 'gex': neg_gex}}
else:
    sorting = {'tot': {'tcr': pos_tcr, 'mhc': pos_mhc, 'hsh': pos_hsh, 'mrk': pos_mrk, 'gex': pos_gex}}

total='tot'
sorting_set = set(list(sorting.keys()) + [total])

# List the keys for barcode files that are custom desgin
custom_barcodes = list()
for brc, relevant in zip(['mhc','hsh','mrk'], [mhc_custom, hsh_custom, mrk_custom]):
    if relevant:
        custom_barcodes.append(brc)

