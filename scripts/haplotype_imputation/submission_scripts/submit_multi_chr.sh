# ARGS1: chromosomes (separated by ,)
# ARGS2: number of threads
# ARGS3: if testing, specify output prefix

chromosomes=$1
genoprefix=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf_multi_chr/preload/no_hypertension
npyprefix=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf_multi_chr/preload/no_hypertension
nthread=$2

outprefix=$3
if [[ -z $outptprefix ]]
then
  outprefix=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf_multi_chr/results/no_hypertension
fi

screen -dmS impute-$i bash -c "bash run_multi_chr.sh $chromosomes $genoprefix $npyprefix $outprefix $nthread"
