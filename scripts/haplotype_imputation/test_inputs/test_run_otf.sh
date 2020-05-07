GENO=/vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr16.h5
SNPyaml=test_inputs/test_snp_list.yaml
PHENO_F=test_inputs/test_father_phenotype.yaml
PHENO_M=test_inputs/test_mother_phenotype.yaml
MODE=basic_em
OUT=test_run_otf.tsv.gz

python impute_po_otf.py \
  --genotype-in-hdf5 $GENO \
  --snp-list-yaml $SNPyaml \
  --father-phenotype-yaml $PHENO_F \
  --mother-phenotype-yaml $PHENO_M \
  --impute-mode $MODE \
  --output $OUT
