
if [[ ! -d logs ]]
then
  mkdir -p logs
fi

for i in `seq 1 22`
do
  screen -dmS ld_clump_$i bash -c "bash submit_ld_clump_one_chromosome.sh \
  <(ls /vol/bmd/yanyul/UKB/neale_lab_gwas/|sed 's#.tsv.bgz##g') \
  $i" 
done

