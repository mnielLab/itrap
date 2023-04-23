#!/usr/bin/env python

dir_fqs = "/home/tuba/herpov/tcr-pmhc-sc-project/data/exp13_RAW/X204SC21092374-Z01-F001/raw_data" # OBS! I have removed the trailing /

#pos_tcr_fq_dir="X204SC21092374-Z01-F001/raw_data/EXP13_TCR1OS"
#neg_tcr_fq_dir="X204SC21092374-Z01-F001/raw_data/EXP13_TCRmix"
#
#pos_mhc_fq_dir="X204SC21092374-Z01-F001/raw_data/EXP13_BC1OS"
#neg_mhc_fq_dir="X204SC21092374-Z01-F001/raw_data/EXP13_BCmix"
#
#pos_hsh_fq_dir="X204SC21092374-Z01-F001/raw_data/EXP13_BC1OS" # Hashing label
#neg_hsh_fq_dir="X204SC21092374-Z01-F001/raw_data/EXP13_BCmix" # Hashing label
#
#pos_mrk_fq_dir="" # Surface/feature marker label
#neg_mrk_fq_dir="" # Surface/feature marker label
#
#pos_gex_fq_dir=""
#neg_gex_fq_dir=""

pos_tcr = "EXP13_TCR1OS-SCI7T010-SCI5T010_HFK5GDSX2"
neg_tcr = "EXP13_TCRmix-SCI7T012-SCI5T012_HFK5GDSX2"

pos_mhc = "EXP13_BC1OS-AK33431-AK33432_HFK5GDSX2"
neg_mhc = "EXP13_BCmix-AK33704-AK33705_HFK5GDSX2"

pos_hsh = "EXP13_BC1OS-AK33431-AK33432_HFK5GDSX2" # Hashing label
neg_hsh = "EXP13_BCmix-AK33704-AK33705_HFK5GDSX2" # Hashing label

pos_mrk = "" # Surface/feature marker label
neg_mrk = "" # Surface/feature marker label

pos_gex = "EXP13_1OS-SCI7T013-SCI5T013_H22JHDSX3"
neg_gex = "EXP13_1OS-SCI7T013-SCI5T013_H22JHDSX3"

mhc_custom = False
hsh_custom = False
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

