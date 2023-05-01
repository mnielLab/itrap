# Improved T cell Receptor Antigen Pairingthrough data-driven filtering of sequencing information from single-cells
This repository contains the code to identify ITRAP filters for single-cell immune profiling data.

## License 
ITRAP is developed by Morten Nielsen's group at the Technical University of Denmark (DTU).
ITRAP code and data can be used freely by academic groups for non-commercial purposes.
If you plan to use ITRAP or any data provided with the script in any for-profit application, you are required to obtain a separate license (contact Morten Nielsen, morni@dtu.dk). 

For scientific questions, please contact Morten Nielsen (mniel@dtu.dk).

## Run ITRAP
Designed with snakemake workflow (v5.7.4)

`snakemake --config exp=exp13 run=run1 --use-conda`

### Run individual steps
Each script may also be run by command line.
For help run scripts/<scriptname> -h
The Snakefile links the required environments for each script.
The environment files are found in envs/.

### Requirements
Anaconda or other Python source (Python 3.7.3)
Specific requirements for each script are logged in envs/

## Data
The pipeline expects a TSV file indexed by 10x barcodes, i.e. GEMs, containing features of TCR, pMHC, & cell hashing.
These features may be generated using [Cellranger multi](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/installation).

The pipeline expects database of TCR-pMHC annotated sequences, which is stored in tools/tcr_dbs.csv.gz.

## Citation


