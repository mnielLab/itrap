import glob
import os
import re
import yaml

"""
Run pipeline from location of the Snakefile:
snakemake -s Snakefile --configfile config/config.13.yaml --config run=1 --cores 4 --use-conda
"""

print(config['run'])

WRK_DIR = workflow.basedir
EXP_DIRECTORY = os.path.join(WRK_DIR, "experiments", config["exp"], config["run"])
RAW_DIRECTORY = os.path.join(WRK_DIR, "data", config["exp"] + "_RAW") #os.path.join(WRK_DIR, "experiments", config["exp"], "raw")
LIB_DIRECTORY = os.path.join(EXP_DIRECTORY, "lib")
MHC_DIRECTORY = os.path.join(EXP_DIRECTORY, "brc")
TCR_DIRECTORY = os.path.join(EXP_DIRECTORY, "tcr")
CAT_DIRECTORY = os.path.join(EXP_DIRECTORY, "cat")
PLT_DIRECTORY = os.path.join(EXP_DIRECTORY, "plt")

include: LIB_DIRECTORY + "/fq_files.py"

##########################
#         Output         #
##########################
TARGET = dict()

include: "workflow/rules/run_cellranger.smk"
if any((mhc_custom, hsh_custom, mrk_custom)):
    include: "workflow/rules/run_brc_pipeline.smk"
include: "workflow/rules/run_cat_pipeline.smk"
include: "workflow/rules/run_analysis.smk"
#include: "workflow/rules/run_pre_plots.smk"
include: "workflow/rules/run_res_plots.smk"

rule all:
    input: TARGET.values(),


