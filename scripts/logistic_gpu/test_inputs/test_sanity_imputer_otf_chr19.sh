HDF5=/vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr19.h5
CHUNKSIZE=30  # 30
NTHREADS=4
PHENO_F=test_inputs/test_phenotype_father.yaml
COVAR_F=test_inputs/test_covariate.yaml
PROBZ=test_inputs/test_prob_z_imputer_otf_chr19.yaml

# PROBZ=test_inputs/test_prob_z_flip.yaml
OUT=test_sanity_imputer_otf_chr19.npy


python run_haplo_logistic_solver.py \
  --genotype-in-hdf5 $HDF5 \
  --variant-chunk-size $CHUNKSIZE \
  --n-threads $NTHREADS \
  --gpu-index 1 \
  --father-phenotype-yaml $PHENO_F \
  --father-covariate-yaml $COVAR_F \
  --haplotype-imputation-yaml $PROBZ \
  --out-npy $OUT
