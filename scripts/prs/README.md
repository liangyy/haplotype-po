Calculate polygenic risk score from GWAS and target genotypes.

Generate YAML for `naive_prs.py` by chromosome.

```
mkdir -p configs
python generate_template.py \
  --yaml-template template.yaml \
  --gwas-list <(cat ../../analysis_output/summary.gwas_neale_lab.csv |tail -n +2 | cut -d, -f 5|awk -F"\/" '{print $6}'| awk -F"?" '{print $1}' | sed 's#.tsv.bgz##g') \
  --chr-num 16 \
  --out-yaml configs/chr16.yaml
```
