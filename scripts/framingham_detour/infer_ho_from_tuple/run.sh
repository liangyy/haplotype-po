# ARGS1: config file 
# ARGS2: log file
# ARGS3: number of cores 

source ~/conda_init.sh
source ~/snmk_init.sh
conda activate haplotype_po

$SNMK --configfile $1 -s infer.snmk --cores $3 -p > $2 2>&1

