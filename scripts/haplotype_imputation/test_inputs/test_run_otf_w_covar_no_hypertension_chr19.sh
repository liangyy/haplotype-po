GENO=/vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr19.h5
SNPyaml=test_inputs/test_full_traits_snp_list_chr19.yaml
SNPcache=test_inputs/test_full_traits_snp_list_chr19.pgz
PHENO_F=test_inputs/test_father_phenotype_no_hypertension.yaml
PHENO_M=test_inputs/test_mother_phenotype_no_hypertension.yaml
COVAR=test_inputs/test_shared_covariate.yaml
MODE=basic_em
OUT=test_run_otf_w_covar_no_hypertension_chr19.tsv.gz
NTHREAD=8

python impute_po_otf.py \
  --genotype-in-hdf5 $GENO \
  --snp-list-yaml $SNPyaml \
  --snp-list-cache $SNPcache \
  --father-phenotype-yaml $PHENO_F \
  --mother-phenotype-yaml $PHENO_M \
  --shared-covariate-yaml $COVAR \
  --impute-mode $MODE \
  --nthread $NTHREAD \
  --output $OUT > test_run_otf_w_covar_no_hypertension_chr19.log 2>&1
