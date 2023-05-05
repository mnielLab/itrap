# Description

This directory contains raw data and data filtered at 3 levels:
- raw: GEMs containing both TCR and pMHC reads. TCR CDR3 sequences must be IUPAC encoded.
- opt_thr: raw data filtered by optimal UMI thresholds defined by ITRAP.
- hla_match: opt_thr data additionally filtered by matching HLA ofpMHC with donor HLA profile.
- tcr: hla_match data additionally filtered to only contain GEMs with both TCRa- and TCRb-chain.


## Contents

.. list-table:: Columns of input data
    :widths: 25 75
    :header-rows: 1
    
    * - Column names
      - Value description
    * - gem
      - Unique identifier for each droplet capture (10x barcode)
    * - clonotype
      - Identifier for each clonotype (annotated by 10x)
    * - ct
      - Identifier for each clonotype (modified from 10x and utilized by ITRAP)
    * - genes_TRA
      - The set of TCRa-chain genes with highest UMI count in a GEM (annotated by 10x)
    * - genes_lst_TRA
      - List of all TCRa-chain genes detected in a GEM (genes sep=";", multiplets sep="|")
    * - genes_TRB
      - The set of TCRb-chain genes with highest UMI count in a GEM (annotated by 10x)
    * - genes_lst_TRB
      - List of all TCRb-chain genes detected in a GEM (genes sep=";", multiplets sep="|")
    * - cdr1_TRA
      - The a-chain CDR1 sequence derived from genes_TRA
    * - cdr2_TRA
      - The a-chain CDR2 sequence derived from genes_TRA
    * - cdr3_TRA
      - The a-chain CDR3 sequence derived from genes_TRA
    * - cdr3_lst_TRA
      - List of all a-chain CDR3 sequences detected in a GEM (sep="|")
    * - cdr1_TRA
      - The a-chain CDR1 sequence derived from genes_TRA
    * - cdr2_TRA
      - The a-chain CDR2 sequence derived from genes_TRA
    * - cdr3_TRB
      - The b-chain CDR3 sequence derived from genes_TRB
    * - cdr3_lst_TRB
      - List of all b-chain CDR3 sequences detected in a GEM (sep="|")
    * - umi_count_TRA
      - UMI count of most abundant TCRa-chain in a GEM
    * - umi_count_lst_TRA
      - List of all TCRa UMI detected in a GEM (sep="|")
    * - delta_umi_TRA
      - Ratio between most abundant UMI and second most abundant UMI, d = UMI_max/(UMI_2+0.25)
    * - umi_count_TRB
      - UMI count of most abundant TCRb-chain in a GEM
    * - umi_count_lst_TRB
      - List of all TCRb UMI detected in a GEM (sep="|")
    * - delta_umi_TRB
      - Ratio between most abundant UMI and second most abundant UMI, d = UMI_max/(UMI_2+0.25)
    * - epitope_rank
      - Rank of peptides (only relevant for plotting)
    * - peptide_HLA
      - peptide and cognate MHC molecule used for staining
    * - HLA_mhc
      - HLA allele from peptide_HLA
    * - peptide_HLA_lst
      - List of all pMHC barcodes detected in a GEM (sep="|")
    * - umi_count_mhc
      - UMI count of most abundant pMHC in a GEM
    * - umi_count_lst_mhc
      - List of all pMCH UMI detected in a GEM (sep="|")
    * - delta_umi_mhc
      - Ratio between most abundant UMI and second most abundant UMI, d = UMI_max/(UMI_2+0.25)
    * - umi_count_mhc_rel
      - Normalized UMI counts
    * - sample_id
      - Identifier for donor given by cell-hashing barcode
    * - sample_id_lst
      - List of all cell-hashing barcodes detected in a GEM (sep="|")
    * - umi_count_cd8
      - UMI count of most abundant cell-hashing barcode in a GEM
    * - umi_count_lst_cd8
      - List of all cell-hashing UMI detected in a GEM (sep="|")
    * - delta_umi_cd8
      - Ratio between most abundant UMI and second most abundant UMI, d = UMI_max/(UMI_2+0.25)
    * - HLA_cd8
      - HLA profile of donor
    * - HLA_lst_cd8
      - List of all HLA alleles associated with the detected sample IDs in a GEM (alleles sep=";", multiplets sep="|")
    * - HLA_match
      - Boolean label for match between HLA of pMHC and donor HLA profile
    * - hto_classification
      - Classification result, with doublets/multiplets named by the top two highest hashtags (Seurat HTO analysis)
    * - hto_global_class
      - Global classification result (singlet, doublet or negative) (Seurat HTO analysis)
    * - hash_id
      - Classification result where doublet IDs are collapsed (Seurat HTO analysis)
    * - hto_max_id
      - Name of hashtag with the highest signal (Seurat HTO analysis)
    * - valid_ct
      - Boolean label for whether clonotype was used for ITRAP threshold optimization
    * - ct_pep
      - Indication of expected binder for given clonotype
    * - gex
      - Boolean label for whether GEM can be concidered to contain a viable cell according to gene expression data