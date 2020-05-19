# ARGS1: genotype in HDF5
# ARGS2: snp list yaml (prefix)
# ARGS3: father pheno yaml
# ARGS4: mother pheno yaml
# ARGS5: covariate yaml
# ARGS6: output prefix
# ARGS7: nthreads

GENO=$1
SNPyaml=$2.yaml
SNPcache=$2.pgz
PHENO_F=$3
PHENO_M=$4
COVAR=$5
MODE=basic_em
OUT=$6.tsv.gz
PKLOUT=$6.pkl.gz
NTHREAD=$7

python impute_po_otf.py \
  --genotype-in-hdf5 $GENO \
  --snp-list-yaml $SNPyaml \
  --snp-list-cache $SNPcache \
  --father-phenotype-yaml $PHENO_F \
  --mother-phenotype-yaml $PHENO_M \
  --shared-covariate-yaml $COVAR \
  --impute-mode $MODE \
  --nthread $NTHREAD \
  --output $OUT \
  --imputer-output $PKLOUT > test_run_otf_per_snp_w_covar_no_hypertension_chr16.log 2>&1
