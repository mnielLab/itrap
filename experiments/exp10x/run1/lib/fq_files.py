#!/usr/bin/env python

dir_fqs = "/home/tuba/herpov/tcr-pmhc-sc-project/experiments/exp10x/raw" # OBS! I have removed the trailing /

pos_tcr = "fake_string"
neg_tcr = ""

pos_mhc = "fake_string"
neg_mhc = ""

pos_hsh = "" # Hashing label
neg_hsh = "" # Hashing label

pos_mrk = "" # Surface/feature marker label
neg_mrk = "" # Surface/feature marker label

pos_gex = ""
neg_gex = ""

mhc_custom = False
hsh_custom = None #False
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

