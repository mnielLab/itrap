import glob
import os
import re
import yaml

#configfile: "config/config.13.yaml" #"config.yaml"
#config['run'] = 'run2'
print(config['run'])

# Load all config variables to be global variables of present script
#f = open("10X.config.yaml")
#globals().update(yaml.load(f)) 
#f.close()
# Or create an object from the content of yaml files
# https://stackoverflow.com/questions/1305532/convert-nested-python-dict-to-object
# https://github.com/Infinidat/munch

WRK_DIR = workflow.basedir
EXP_DIRECTORY = os.path.join(WRK_DIR, "experiments", config["exp"], config["run"])
RAW_DIRECTORY = os.path.join(WRK_DIR, "data", config["exp"] + "_RAW") #os.path.join(WRK_DIR, "experiments", config["exp"], "raw")
LIB_DIRECTORY = os.path.join(EXP_DIRECTORY, "lib")
MHC_DIRECTORY = os.path.join(EXP_DIRECTORY, "brc") # change to BRC
TCR_DIRECTORY = os.path.join(EXP_DIRECTORY, "tcr")
CAT_DIRECTORY = os.path.join(EXP_DIRECTORY, "cat")
PLT_DIRECTORY = os.path.join(EXP_DIRECTORY, "plt")



# include or load as yaml/dict?
# https://github.com/tanaes/snakemake_assemble/blob/master/bin/snakefiles/folders
#with open(LIB_DIRECTORY + "/path_lib.yaml", 'r') as stream:
#    path = yaml.safe_load(stream)

include: LIB_DIRECTORY + "/fq_files.py"

##########################
#         Output         #
##########################
TARGET = dict()

#TARGET['readme'] = EXP_DIRECTORY + "/README.md"

include: "workflow/rules/run_cellranger.smk"
if any((mhc_custom, hsh_custom, mrk_custom)):
    include: "workflow/rules/run_brc_pipeline.smk"
include: "workflow/rules/run_cat_pipeline.smk"
include: "workflow/rules/run_analysis.smk"
#include: "workflow/rules/run_pre_plots.smk"
include: "workflow/rules/run_res_plots.smk"

rule all:
    input: TARGET.values(), #README, RANGER, CONTIG, TESTER


