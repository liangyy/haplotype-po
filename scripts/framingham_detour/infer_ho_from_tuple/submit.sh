# ARGS1: which rule in snmk to apply. If not apply it will follow the default behavior of snakemake

screen -dmS from_pedigree bash -c "bash run.sh config.yaml run.log 6 $1"

