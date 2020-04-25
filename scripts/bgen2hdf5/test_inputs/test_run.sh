DATADIR=test_inputs
CHR=16

GWASYAML=$DATADIR/test_gwas.yaml
SNPMAP=/vol/bmd/yanyul/UKB/haplotype_imputation/snp_map_for_neale_lab_gwas.with_sign.tsv.gz
BGEN=/vol/bmd/meliao/data/haplotype/hap/ukb_hap_chr$CHR\_v2.bgen
BGI=/vol/bmd/meliao/data/haplotype/hap_bgi/ukb_hap_chr$CHR\_v2.bgen.bgi
SAMPLE=/vol/bmd/meliao/data/haplotype/link_files/ukb1952_v2_s487398.sample


python naive_prs.py \
  --bgen $BGEN \
  --bgi $BGI \
  --sample $SAMPLE \
  --output-hdf5 $DATADIR/test_out.h5 \
  --snp-chunk-size 5 \
  --bgen-writing-cache-size 100 \
  --max-sample-chunk-size 10000 \
  --max-snp-chunk-size 4 \
  --first-n-snp 27
  
