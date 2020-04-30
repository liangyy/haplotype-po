HDF5=/vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr16.h5
CHUNKSIZE=100
NTHREADS=4
PHENO=test_inputs/test_phenotype_father.yaml
COVAR=test_inputs/test_covariate.yaml
INDIVLIST=tmp_sanity_sum_indiv_list.txt
cat ../../analysis_output/test_em_output.tsv|cut -f 1|tail -n +2 > $INDIVLIST
OUT=test_sanity_sum.npy


python run_logistic_solver.py \
  --genotype-in-hdf5 $HDF5 \
  --variant-chunk-size $CHUNKSIZE \
  --n-threads $NTHREADS \
  --gpu-index 0 \
  --phenotype-yaml $PHENO \
  --covariate-yaml $COVAR \
  --individual-list $INDIVLIST \
  --out-npy $OUT
