CHR=22

GWASYAML=test_gwas.yaml
SNPMAP=/vol/bmd/yanyul/UKB/haplotype_imputation/snp_map_for_neale_lab_gwas.tsv.gz
BGEN=/vol/bmd/meliao/data/haplotype/hap/ukb_hap_chr$CHR\_v2.bgen
BGI=/vol/bmd/meliao/data/haplotype/hap_bgi/ukb_hap_chr$CHR\_v2.bgen.bgi
SAMPLE=/vol/bmd/meliao/data/haplotype/link_files/ukb1952_v2_s487398.sample


python naive_prs.py \
  --gwas-yaml $GWASYAML \
  --snp-map $SNPMAP \
  --bgen $BGEN \
  --bgi $BGI \
  --sample $SAMPLE \
  --pval-cutoffs 0.01,0.5,1 \
  --chromosome $CHR \
  --output-hdf5 test_out.h5 
  