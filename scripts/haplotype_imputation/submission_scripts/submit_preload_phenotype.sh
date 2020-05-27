genodir=/lambda_stor/data/yanyul/UKB/ukb_hap_v2_to_hdf5_genomewide
phenof=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/haplotype_imputation/submission_scripts/father_phenotype_no_hypertension.yaml
phenom=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/haplotype_imputation/submission_scripts/mother_phenotype_no_hypertension.yaml 
covar=/lambda_stor/data/yanyul/GitHub/haplotype-po/scripts/haplotype_imputation/submission_scripts/covariates.yaml

outdir=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf_multi_chr

mkdir -p $outdir/preload

screen -dmS preload-pheno bash -c "bash preload_phenotype.sh $genodir/ukb_hap_v2_to_hdf5.chr1.h5 $phenof $phenom $covar $outdir/preload/no_hypertension"