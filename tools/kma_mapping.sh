ARGS=$3 #"-1t1" $(echo -e "${ARGS}" | tr -d '[:space:]')

BARCODE_TEMPLATES=$1 #"/home/tuba/herpov/tcr-pmhc-sc-project/data/exp3_MHC/barcode_library/barcode_templates.fa"
READS=$2 #"/home/tuba/herpov/tcr-pmhc-sc-project/data/exp3_MHC/processed/longranger_clean/IONTORRENT.R1.gems.no_umi.no_adapters.fq"

KMA=/home/tuba-nobackup/shared/src/kma/kma

KMA_DIR=$(echo $4 | rev | cut -d"/" -f3- | rev) 

#mkdir -p $KMA_DIR/output

$KMA index -i $BARCODE_TEMPLATES -o $KMA_DIR/templates
$KMA -i $READS -o $KMA_DIR/output/mapping -t_db $KMA_DIR/templates $ARGS
