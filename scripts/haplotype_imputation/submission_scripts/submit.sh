# ARGS1: from chr number
# ARGS2: to chr number
# ARGS3: number of threads

genodir=/lambda_stor/data/yanyul/UKB/ukb_hap_v2_to_hdf5_genomewide
snpdir=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/haplotype_imputation/submission_scripts/full_traits_snp_list
phenof=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/haplotype_imputation/submission_scripts/father_phenotype_no_hypertension.yaml
phenom=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/haplotype_imputation/submission_scripts/mother_phenotype_no_hypertension.yaml 
covar=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/haplotype_imputation/submission_scripts/covariates.yaml
nthread=$3

outdir=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf


for i in `seq $1 $2`
do
  screen -dmS impute-$i bash -c "bash run.sh $genodir/ukb_hap_v2_to_hdf5.chr$i.h5 $snpdir/chr$i $phenof $phenom $covar $outdir/chr$i $nthread"
done
