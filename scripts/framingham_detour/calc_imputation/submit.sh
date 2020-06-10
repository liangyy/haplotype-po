# ARGS1: number of cores in snakemake call
# ARGS2: rule name to apply (if not specify, will follow snakemake default)

screen -dmS impute_en bash -c "bash run.sh $1 en $2"
screen -dmS impute_dapgw bash -c "bash run.sh $1 dapgw $2"

