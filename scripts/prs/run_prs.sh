# ARGS1: chromosome, e.g. 16
# ARGS2: GWAS YAML
# ARGS3: pval cutoffs
# ARGS4: outdir
# ARGS5: name tag
# ARGS6: gwas list

CHR=$1
OUTDIR=$4
NAMETAG=$5

GWASYAML=$2
PVALS=$3  # "5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1"
GWASLIST=$6
SNPMAP=/vol/bmd/yanyul/UKB/haplotype_imputation/snp_map_for_neale_lab_gwas.with_sign.tsv.gz
BGEN=/vol/bmd/meliao/data/haplotype/hap/ukb_hap_chr$CHR\_v2.bgen
BGI=/vol/bmd/meliao/data/haplotype/hap_bgi/ukb_hap_chr$CHR\_v2.bgen.bgi
SAMPLE=/vol/bmd/meliao/data/haplotype/link_files/ukb1952_v2_s487398.sample

mkdir -p $OUTDIR

python naive_prs.py \
  --gwas-yaml $GWASYAML \
  --snp-map $SNPMAP \
  --bgen $BGEN \
  --bgi $BGI \
  --sample $SAMPLE \
  --pval-cutoffs $PVALS \
  --chromosome $CHR \
  --output-hdf5 $OUTDIR/prs_naive.$NAMETAG.chr$CHR.h5 \
  --snp-chunk-size 100 \
  --cache-path-for-gwas-dict $OUTDIR/prs_naive.chr$CHR.pgz \
  --gwas-list $GWASLIST \
  > $OUTDIR/prs_naive.$NAMETAG.chr$CHR.log 2>&1

