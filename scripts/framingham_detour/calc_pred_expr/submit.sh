# ARGS1: configfile

source ~/conda_init.sh
source ~/snmk_init.sh
conda activate metaxcan

cmd="$SNMK --configfile $1 -s run.snmk -np"  

