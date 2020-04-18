# ARGS: outdir
cat ../../analysis_output/summary.gwas_neale_lab.csv |cut -f 5 -d,|sed 's#"##g' > download_neale_tmp.sh
mydir=`pwd`
mkdir -p $1
cd $1
bash $mydir/download_neale_tmp.sh
rm $mydir/download_neale_tmp.sh

