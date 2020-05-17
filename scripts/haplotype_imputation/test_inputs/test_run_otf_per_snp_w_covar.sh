GENO=/vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr16.h5
SNPyaml=test_inputs/test_snp_list.yaml
SNPcache=test_inputs/test_snp_list.pgz
COVAR=test_inputs/test_shared_covariate.yaml
PHENO_F=test_inputs/test_father_phenotype.yaml
PHENO_M=test_inputs/test_mother_phenotype.yaml
MODE=per_snp_em
OUT=test_run_otf_per_snp_w_covar.tsv.gz

python impute_po_otf.py \
  --genotype-in-hdf5 $GENO \
  --snp-list-yaml $SNPyaml \
  --snp-list-cache $SNPcache \
  --father-phenotype-yaml $PHENO_F \
  --mother-phenotype-yaml $PHENO_M \
  --shared-covariate-yaml $COVAR \
  --impute-mode $MODE \
  --output $OUT
