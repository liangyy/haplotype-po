GENO=/vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr16.h5
SNPyaml=test_inputs/test_full_traits_snp_list.yaml
SNPcache=test_inputs/test_full_traits_snp_list.pgz
PHENO_F=test_inputs/test_father_phenotype.yaml
PHENO_M=test_inputs/test_mother_phenotype.yaml
COVAR=test_inputs/test_shared_covariate.yaml
MODE=basic_em
OUT=test_run_otf_w_covar_full_traits_chr16.tsv.gz

python impute_po_otf.py \
  --genotype-in-hdf5 $GENO \
  --snp-list-yaml $SNPyaml \
  --snp-list-cache $SNPcache \
  --father-phenotype-yaml $PHENO_F \
  --mother-phenotype-yaml $PHENO_M \
  --shared-covariate-yaml $COVAR \
  --impute-mode $MODE \
  --output $OUT > test_run_otf_w_covar_full_traits_chr16.log 2>&1
