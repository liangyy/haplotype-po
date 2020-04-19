GENOpattern=/vol/bmd/meliao/data/haplotype/hap/ukb_hap_chr{chr}_v2.bgen
OUTdir=/vol/bmd/yanyul/UKB/ukb_hap_bgen_reader_metafile/ukb_hap_v2

python build_bgen_reader_meta_file.py \
  --bgen-pattern $GENOpattern \
  --out-prefix $OUTdir
  