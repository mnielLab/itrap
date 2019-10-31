#!/bin/bash



# get illumna seq data for our 10x barcoding pilot study, which was sequenced at DTU bioengineering

# install bs and bs-cp (use curl not wget, i think tubas wget it out of date)


cd /home/people/marquard/tuba/tcr_seq/10x-barcoding/2019_02_25_Exp2/data
mkdir fastq
mkdir raw

bs list datasets

# +-------------+-------------------------------------+---------------------+---------------------+
# |    Name     |                 Id                  |    Project.Name     |   DataSetType.Id    |
# +-------------+-------------------------------------+---------------------+---------------------+
# | TCR_2_L001  | ds.9ab0b98f8fb248b2a74f05c6289c3213 | MS_2x150_TCR_190222 | illumina.fastq.v1.8 |
# | TCR_L001    | ds.6717ab9d40fc414599afbce2d590c8ca | MS_2x150_TCR_181130 | illumina.fastq.v1.8 |
# | input2_L001 | ds.2f8096f0efc545a6a06add4b435fa12b | MS_1x200_MHC_181129 | illumina.fastq.v1.8 |
# | input3_L001 | ds.d1d90581ec2d4dc4935eb8974ce03a7b | MS_1x200_MHC_181129 | illumina.fastq.v1.8 |
# | input1_L001 | ds.db4ff511ce5c4549b02730d9ccc91ad4 | MS_1x200_MHC_181129 | illumina.fastq.v1.8 |
# | MHC_L001    | ds.346d77eab8134dbfbdba0cb12e3701ad | MS_1x200_MHC_181129 | illumina.fastq.v1.8 |
# | CD8_L001    | ds.a6603227f1b94e19bd19c330ea2a6a68 | MS_1x200_MHC_181129 | illumina.fastq.v1.8 |
# +-------------+-------------------------------------+---------------------+---------------------+


# Get RAW intensities to make fastqs ourselves
# | TCR_2_L001  | ds.9ab0b98f8fb248b2a74f05c6289c3213 | MS_2x150_TCR_190222 | illumina.fastq.v1.8 |
bs-cp -v https://basespace.illumina.com/Run/MS_2x150_TCR_190222 .

# Next for barcodes:

# Get FASTQS
# bs download dataset -i ds.2f8096f0efc545a6a06add4b435fa12b -o fastq











# [marquard@tuba 2019_02_25_Exp2]$ bs list project
# +---------------------+-----------+------------+
# |        Name         |    Id     | TotalSize  |
# +---------------------+-----------+------------+
# | MS_1x200_MHC_190227 | 121051942 | 0          |
# | MS_2x150_TCR_190222 | 119093210 | 194345     |
# | MS_2x150_TCR_181130 | 106671604 | 34081      |
# | MS_1x200_MHC_181129 | 106663563 | 2894827534 |
# +---------------------+-----------+------------+



# [marquard@tuba 2019_02_25_Exp2]$ bs list run
# +------------------------------------+-----------+---------------------+----------+
# |                Name                |    Id     |   ExperimentName    |  Status  |
# +------------------------------------+-----------+---------------------+----------+
# | 190227_M02023_0490_000000000-C8825 | 164218054 | MS_1x200_MHC_190227 | Complete |
# | 190222_M02023_0488_000000000-BWKRT | 162759601 | MS_2x150_TCR_190222 | Complete |
# | 181130_M02023_0480_000000000-BW5CW | 144356214 | MS_2x150_TCR_181130 | Complete |
# | 181129_M02023_0479_000000000-C7JYR | 144095952 | MS_1x200_MHC_181129 | Complete |
# +------------------------------------+-----------+---------------------+----------+