# ARGS1: chromosomes (separated by ,)
# ARGS2: genotype prefix pattern
# ARGS3: prefix of other npy preloaded files
# ARGS4: output prefix

CHROMS=$1
genotypePREFIX=$2
preloadNPY=$3
OUT=$4

cd ../

python impute_otf_multi_chr.py \
  --genotype-prefix-pattern $genotypePREFIX.chr{chr_num} \
  --chromosomes $CHROMS \
  --npy-prefix $preloadNPY \
  --output-prefix $OUT.chr \
  --imputer-output $OUT.pkl.gz 
  
