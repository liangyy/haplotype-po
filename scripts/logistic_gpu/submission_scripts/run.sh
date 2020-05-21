# ARGS1: input hdf5 genotype
# ARGS2: phenotype of father (this script only do one parent)
# ARGS3: covariate of father
# ARGS4: imputation table
# ARGS5: output prefix
# ARGS6: GPU index

# fixed meta parameters
CHUNKSIZE=30  
NTHREADS=4

# inputs
HDF5=$1
PHENO_F=$2
COVAR_F=$3
PROBZ=$4
GPUindex=$6

# generate prob z YAML
PROBZyaml=$5.prob_z.yaml
cat prob_z_template.yaml |sed "s#PLACEHOLDER#$PROBZ#" > $PROBZyaml

# output
OUT=$5.npy
MYLOG=$5.log

cd ../

python run_haplo_logistic_solver.py \
  --genotype-in-hdf5 $HDF5 \
  --variant-chunk-size $CHUNKSIZE \
  --n-threads $NTHREADS \
  --gpu-index $GPUindex \
  --father-phenotype-yaml $PHENO_F \
  --father-covariate-yaml $COVAR_F \
  --haplotype-imputation-yaml $PROBZyaml \
  --out-npy $OUT > $MYLOG 2>&1
