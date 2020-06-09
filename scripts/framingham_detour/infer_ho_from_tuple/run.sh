# ARGS1: config file 
# ARGS2: log file
# ARGS3: number of cores 
# ARGS4: which rule in snmk to apply (optional). If not apply, will follow default behavior of snakemake

source ~/conda_init.sh
source ~/snmk_init.sh
conda activate haplotype_po

$SNMK --configfile $1 -s infer.snmk $4 --cores $3 -p > $2 2>&1

