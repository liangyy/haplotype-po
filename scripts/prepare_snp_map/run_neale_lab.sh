# inputs
GENOpattern=/vol/bmd/meliao/data/haplotype/hap/ukb_hap_chr{chr}_v2.bgen
SAMPLEfile=/vol/bmd/meliao/data/haplotype/link_files/ukb1952_v2_s487398.sample
GWASfile=/vol/bmd/yanyul/UKB/haplotype_imputation/20002_1065.gwas.imputed_v3.both_sexes.tsv.bgz

# output
OUTfile=/vol/bmd/yanyul/UKB/haplotype_imputation/snp_map_for_neale_lab_gwas.tsv.gz

# run
python build_snp_map_for_neale_lab_gwas.py \
  --genotype-pattern $GENOpattern \
  --genotype-sample $SAMPLEfile \
  --output $OUTfile \
  --gwas $GWASfile > \
  run_neale_lab.log 2>&1
  