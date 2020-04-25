for GRP in `seq 1 10`
do
  screen -dmS PRS$GRP bash -c "run_prs.sh 16 /vol/bmd/yanyul/GitHub/haplotype-po/scripts/prs/configs/chr16.yaml 5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1 $outdir group$GRP"
done