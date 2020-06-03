# ARGS1: configfile
# ARGS2: jobname

source ~/conda_init.sh
source ~/snmk_init.sh
conda activate metaxcan

cmd="$SNMK --configfile $1 -s run.snmk -p > $2.log 2>&1"  
screen -dmS $2 bash -c "$cmd"
