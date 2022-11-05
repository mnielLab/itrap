import yaml
import os

rule construct_template_database:
	input:
		brc_info = os.path.join(LIB, 'barcode_design.yml'),
		nomencl = os.path(join(LIB, 'barcode_specificity_annotations.xlsx'))
	output:
		MHC_DIRECTORY + "/barcode_library/barcode_templates.fa"
	conda:
		"snake_envs/biopython.yaml"
	script:
		"./scripts/A2_construct_template_database.py"
