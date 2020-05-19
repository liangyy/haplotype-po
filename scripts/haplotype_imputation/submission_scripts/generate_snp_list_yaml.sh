# generate snp list yaml

mkdir -p full_traits_snp_list

for i in `seq 1 22`
do
  python ../generate_snp_list.py --yaml-template <(cat template.yaml | sed "s/XXX/$i/g") --trait-yaml full_traits.yaml --out-yaml full_traits_snp_list/chr$i.yaml
done