Calculate polygenic risk score from GWAS and target genotypes.

Generate YAML for `naive_prs.py` by chromosome.

```
mkdir -p configs
python generate_yaml.py \
  --yaml-template template.yaml \
  --gwas-list <(cat ../../analysis_output/summary.gwas_neale_lab.csv |tail -n +2 | cut -d, -f 5|awk -F"\/" '{print $6}'| awk -F"?" '{print $1}' | sed 's#.tsv.bgz##g') \
  --chr-num 16 \
  --out-yaml configs/chr16.yaml
```

Run `naive_prs.py`

```
outdir=/vol/bmd/yanyul/UKB/haplotype_imputation/prs
bash run_prs.sh \
  16 \
  /vol/bmd/yanyul/GitHub/haplotype-po/scripts/prs/configs/chr16.yaml \
  5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,0.01,0.05,0.1,0.5,1 \
  $outdir \
  name-of-the-run \
  <(cat ../../analysis_output/analysis_batch.gwas_neale_lab.txt |grep _group1_ | cut -f 1)
```
