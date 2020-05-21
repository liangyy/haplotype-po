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
MYLOG=$6.log
NTHREAD=$7

cd ../

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
  --imputer-output $PKLOUT \
  --debug-cache-prefix $6 > $MYLOG 2>&1

