#!/usr/bin/env python

dir_fqs = "/home/tuba/herpov/tcr-pmhc-sc-project/data/exp13_RAW/X204SC21092374-Z01-F001/raw_data" # OBS! I have removed the trailing /

pos_tcr = "EXP13_TCR2OS-SCI7T011-SCI5T011_HFK5GDSX2"
neg_tcr = "EXP13_TCRmix-SCI7T012-SCI5T012_HFK5GDSX2"

pos_mhc = "EXP13_BC2OS_costum-1_HFK5GDSX2" #"fake_reads" #
neg_mhc = "EXP13_BCmix-AK33704-AK33705_HFK5GDSX2"

pos_hsh = "EXP13_BC2OS-AK33706-AK33707_HFK5GDSX2" # Hashing label
neg_hsh = "EXP13_BCmix-AK33704-AK33705_HFK5GDSX2" # Hashing label

pos_mrk = "" # Surface/feature marker label
neg_mrk = "" # Surface/feature marker label

pos_gex = "EXP13_2OS-SCI7T014-SCI5T014_H22JHDSX3"
neg_gex = ""

mhc_custom = True
hsh_custom = False
mrk_custom = None

#sorting = {'pos': {'tcr': pos_tcr, 'mhc': pos_mhc, 'hsh': pos_hsh, 'mrk': pos_mrk, 'gex': pos_gex},
#           'neg': {'tcr': neg_tcr, 'mhc': neg_mhc, 'hsh': neg_hsh, 'mrk': neg_mrk, 'gex': neg_gex}}

if any((neg_tcr, neg_mhc, neg_hsh, neg_mrk, neg_gex)):
    sorting = {'pos': {'tcr': pos_tcr, 'mhc': pos_mhc, 'hsh': pos_hsh, 'mrk': pos_mrk, 'gex': pos_gex},
               'neg': {'tcr': neg_tcr, 'mhc': neg_mhc, 'hsh': neg_hsh, 'mrk': neg_mrk, 'gex': neg_gex}}
else:
    sorting = {'tot': {'tcr': pos_tcr, 'mhc': pos_mhc, 'hsh': pos_hsh, 'mrk': pos_mrk, 'gex': pos_gex}}

total = 'tot'
sorting_set = set(list(sorting.keys()) + [total])

# List the keys for barcode files that are custom desgin
custom_barcodes = list()
for brc, relevant in zip(['mhc','hsh','mrk'], [mhc_custom, hsh_custom, mrk_custom]):
    if relevant:
        custom_barcodes.append(brc)


