PRS=/vol/bmd/yanyul/UKB/haplotype_imputation/prs/prs_naive.full.chr16.h5
PRSyaml=test_inputs/test_map_prs.yaml
PHENO_F=test_inputs/test_father_phenotype.yaml
PHENO_M=test_inputs/test_mother_phenotype.yaml
MODE=basic_em
OUT=test_run.tsv.gz

python impute_parent_of_origin.py \
  --prs-matrix $PRS \
  --map-prs-yaml $PRSyaml \
  --father-phenotype-yaml $PHENO_F \
  --mother-phenotype-yaml $PHENO_M \
  --impute-mode $MODE \
  --output $OUT
