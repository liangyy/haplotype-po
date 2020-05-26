GENO=/vol/bmd/yanyul/UKB/ukb_hap_v2_to_hdf5/ukb_hap_v2_to_hdf5.chr16.h5
SNPyaml=test_inputs/test_snp_list.yaml
SNPcache=test_inputs/test_snp_list.pgz
COVAR=test_inputs/test_shared_covariate.yaml
PHENO_F=test_inputs/test_preload_father_phenotype.yaml
PHENO_M=test_inputs/test_preload_mother_phenotype.yaml
OUTPRE=test_preload_otf

python impute_otf_preload.py \
  --genotype-in-hdf5 $GENO \
  --father-phenotype-yaml $PHENO_F \
  --mother-phenotype-yaml $PHENO_M \
  --shared-covariate-yaml $COVAR \
  --load-mode phenotype \
  --output-prefix $OUTPRE 

python impute_otf_preload.py \
  --genotype-in-hdf5 $GENO \
  --snp-list-yaml $SNPyaml \
  --snp-list-cache $SNPcache \
  --indiv-npy $OUTPRE.individual_id.npy \
  --pheno-npy $OUTPRE.phenotype.npy \
  --load-mode genotype \
  --output-prefix $OUTPRE 

