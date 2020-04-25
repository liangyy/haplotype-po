# ARGS1: chromosome, e.g. 16
# ARGS2: outdir

CHR=$1
OUTDIR=$2

BGEN=/vol/bmd/meliao/data/haplotype/hap/ukb_hap_chr$CHR\_v2.bgen
BGI=/vol/bmd/meliao/data/haplotype/hap_bgi/ukb_hap_chr$CHR\_v2.bgen.bgi
SAMPLE=/vol/bmd/meliao/data/haplotype/link_files/ukb1952_v2_s487398.sample

mkdir -p $OUTDIR

python run_ukb_hap_bgen_to_hdf5.py \
  --bgen $BGEN \
  --bgi $BGI \
  --sample $SAMPLE \
  --output-hdf5 $OUTDIR/ukb_hap_v2_to_hdf5.chr$CHR.h5 \
  --snp-chunk-size 100 \
  --bgen-writing-cache-size 1000 \
  > $OUTDIR/ukb_hap_v2_to_hdf5.chr$CHR.log 2>&1

