# ARGS1: nthreads

source ~/conda_init.sh
source ~/snmk_init.sh
conda activate haplotype_po

wkdir=`pwd`

cd ../
MYLOG=$wkdir/submit.log

screen -dmS per_chr_ldsc bash -c "$SNMK --configfile $wkdir/config.yaml -s per_chr_h2.snmk --cores $1 -p > $MYLOG 2>&1"

