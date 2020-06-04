# ARGS1: config file
# ARGS2: pipeline dir
# ARGS3: log file
# ARGS4: number of cores

pipedir=$2
myconfig=$1

source ~/conda_init.sh
source ~/snmk_init.sh
conda activate haplotype_po

wkdir=`pwd`
cd $pipedir

$SNMK -s pca.snmk --configfile $wkdir/$myconfig -p > $wkdir/$3 --cores $4 2>&1


