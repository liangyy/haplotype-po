# ARGS1: config file
# ARGS2: pipeline dir
# ARGS3: log file

pipedir=$2
myconfig=$1

source ~/conda_init.sh
source ~/snmk_init.sh
conda activate haplotype_po

wkdir=`pwd`
cd $pipedir

$SNMK -s peer.snmk --configfile $wkdir/$myconfig -p > $wkdir/$3 2>&1


