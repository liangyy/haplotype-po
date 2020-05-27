INPREFIX=test_preload_otf

OUT=test_multi_chr16_only.tsv.gz
OUTPKL=test_multi_chr16_only.pkl.gz

python impute_po_otf.py \
  --genotype-prefix-pattern $INPREFIX.chr{chr_num} \
  --chromosomes 16 \
  --npy-prefix $INPREFIX \
  --output $OUT \
  --imputer-output $OUTPKL 
  