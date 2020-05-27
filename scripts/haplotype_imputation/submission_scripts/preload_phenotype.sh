# ARGS1: genotype in HDF5
# ARGS2: father pheno yaml
# ARGS3: mother pheno yaml
# ARGS4: covariate yaml
# ARGS5: output prefix

# To load the individual list
GENO=$1

PHENO_F=$2
PHENO_M=$3
COVAR=$4
OUTPRE=$5

cd ../

python impute_otf_preload.py \
  --genotype-in-hdf5 $GENO \
  --father-phenotype-yaml $PHENO_F \
  --mother-phenotype-yaml $PHENO_M \
  --shared-covariate-yaml $COVAR \
  --load-mode phenotype \
  --output-prefix $OUTPRE 
  