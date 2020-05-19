# ARGS1: from chr number
# ARGS2: to chr number
# ARGS3: number of threads

genodir=/lambda_stor/data/yanyul/UKB/ukb_hap_v2_to_hdf5_genomewide
snpdir=full_traits_snp_list
phenof=father_phenotype_no_hypertension.yaml
phenom=mother_phenotype_no_hypertension.yaml 
covar=covariates.yaml
nthread=$3

outdir=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf


for i in `seq $1 $2`
do
  screen -dmS impute-$i "bash -c "run.sh $genodir/ukb_hap_v2_to_hdf5.chr$i.h5 $snpdir/chr$i $phenof $phenom $covar $/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf/chr$i $nthread
done