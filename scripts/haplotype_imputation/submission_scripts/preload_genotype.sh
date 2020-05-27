# ARGS1: genotype in HDF5
# ARGS2: snp list yaml (prefix)
# ARGS3: prefix of corresponding preloaded phenotypes 
# ARGS4: output prefix

# To load the individual list
GENO=$1
SNPyaml=$2.yaml
SNPcache=$2.pgz
INPRE=$3
OUTPRE=$4

cd ../

python impute_otf_preload.py \
  --genotype-in-hdf5 $GENO \
  --snp-list-yaml $SNPyaml \
  --snp-list-cache $SNPcache \
  --indiv-npy $INPRE.individual_id.npy \
  --pheno-npy $INPRE.phenotype.npy \
  --load-mode genotype \
  --output-prefix $OUTPRE > $OUTPRE.log 2>&1
  
