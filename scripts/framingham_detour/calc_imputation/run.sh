# ARGS1: number of cores
# ARGS2: name_tag 
# ARGS3: rule name

source ~/conda_init.sh
source ~/snmk_init.sh
conda activate haplotype_po


$SNMK --configfile config.yaml -s impute.snmk -p $3 --config name_tag=$2 --cores $1 > run_$2.log 2>&1


