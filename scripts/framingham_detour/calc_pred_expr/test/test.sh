conda activate metaxcan

SCRIPTPATH=/lambda_stor/data/yanyul/GitHub/MetaXcan/software/Predict.py

predDB=/lambda_stor/data/yanyul/misc_data/gtex_v8_pred_models/eqtl/elastic_net_models_dapgw/dapgw_Whole_Blood.db
myVCF=/lambda_stor/data/yanyul/Framingham/imputed_hrc1.1/chr22.dose.vcf.gz
LOchain=/lambda_stor/data/yanyul/misc_data/liftover_chains/hg19ToHg38.over.chain.gz

python $SCRIPTPATH \
  --model_db_path $predDB \
  --model_db_snp_key varID \
  --vcf_genotypes $myVCF \
  --vcf_mode haplotyped \
  --liftover $LOchain \
  --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
  --prediction_output test_out.haplotyped.txt \
  --prediction_summary_output test_summary.haplotyped.txt \
  --verbosity 9 \
  --throw 

python $SCRIPTPATH \
  --model_db_path $predDB \
  --model_db_snp_key varID \
  --vcf_genotypes $myVCF \
  --vcf_mode genotyped \
  --liftover $LOchain \
  --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
  --prediction_output test_out.genotyped.txt \
  --prediction_summary_output test_summary.genotyped.txt \
  --verbosity 9 \
  --throw

conda deactivate
conda activate haplotype_po
jupyter nbconvert --to notebook  --inplace --execute test_check_result.ipynb
