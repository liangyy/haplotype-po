HDF5=/vol/bmd/yanyul/GitHub/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr16.h5
CHUNKSIZE=100
NTHREADS=4
PHENO=test_inputs/test_phenotype.yaml
COVAR=test_inputs/test_covariate.yaml
INDIVLIST='zcat ukb_eur.father.csv.gzip | grep -v "\-2"|cut -f 1 -d,|tail -n +2'
OUT=test_run.npy

python run_logistic_solver.py \
  --genotype-in-hdf5 $HDF5 \
  --variant-chunk-size $CHUNKSIZE \
  --n-threads $NTHREADS \
  --gpu-index 1 \
  --phenotype-yaml $PHENO \
  --covariate-yaml $COVAR \
  --individual-list <($INDIVLIST) \
  --out-npy $OUT