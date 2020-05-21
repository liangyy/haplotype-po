# ARGS1: list of chr
# ARGS2: list of gpu index
# ARGS3: prefix of imputation files
# ARGS4: suffix of imputation files

# args
IFS=',' read -r -a chrArray <<< $1
IFS=',' read -r -a gpuArray <<< $2


# some inputs
genodir=/lambda_stor/data/yanyul/UKB/ukb_hap_v2_to_hdf5_genomewide
phenof=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/logistic_gpu/submission_scripts/father_phenotypes.yaml
covar=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/logistic_gpu/submission_scripts/covariates.yaml
impPrefix=$3
impSuffix=$4

outdir=/lambda_stor/data/yanyul/UKB/haplotype_imputation/gwas

for index in "${!chrArray[@]}"
do
  chrom="${chrArray[index]}"
  gpuidx="${gpuArray[index]}"
  screen -dmS gwas-$chrom bash -c "bash run.sh $genodir/ukb_hap_v2_to_hdf5.chr$chrom.h5 $phenof $covar $impPrefix$chrom$impSuffix $outdir/chr$chrom $gpuidx"
done

