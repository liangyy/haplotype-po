cd /vol/bmd/yanyul/GitHub/haplotype-po/scripts/bgen2hdf5/ 
for chr in `seq 1 22`
do
  screen -dmS convert_$chr bash -c "bash run_bgen2hdf5.sh $chr /vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5_genomewide"
done
