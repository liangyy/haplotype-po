# ARGS1: from chr number
# ARGS2: to chr number
# ARGS3: preloaded phenotype prefix

genodir=/lambda_stor/data/yanyul/UKB/ukb_hap_v2_to_hdf5_genomewide
snpdir=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/haplotype_imputation/submission_scripts/no_hypertension_snp_list
inprefix=$3

outdir=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf_multi_chr/preload


for i in `seq $1 $2`
do
  screen -dmS preload-$i bash -c "bash preload_genotype.sh $genodir/ukb_hap_v2_to_hdf5.chr$i.h5 $snpdir/chr$i $inprefix $outdir/no_hypertension.chr$i"
done

  
