INPREFIX=test_preload_otf

OUT=test_multi_chr16_only
OUTPKL=test_multi_chr16_only.pkl.gz

python impute_otf_multi_chr.py \
  --genotype-prefix-pattern $INPREFIX.chr{chr_num} \
  --chromosomes 16 \
  --npy-prefix $INPREFIX \
  --output-prefix $OUT \
  --imputer-output $OUTPKL 
  
