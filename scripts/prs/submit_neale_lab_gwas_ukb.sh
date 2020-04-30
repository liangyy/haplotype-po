# ARGS1: outdir

outdir=$1
screen -dmS PRS_full bash -c "bash run_prs.sh 16 /vol/bmd/yanyul/GitHub/haplotype-po/scripts/prs/configs/chr16.yaml 5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1 $outdir full"


