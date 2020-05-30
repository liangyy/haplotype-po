inprefix=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf_multi_chr/results/no_hypertension
insuffix=.tsv.gz
outprefix=/lambda_stor/data/yanyul/UKB/haplotype_imputation/haplotype_impute_otf_multi_chr/results/no_hypertension
outsuffix=.w_baseline.tsv.gz

for i in `seq 1 22`
do
  python add_baseline.py --input $inprefix$i$insuffix --output $outprefix$i$outsuffix
done
