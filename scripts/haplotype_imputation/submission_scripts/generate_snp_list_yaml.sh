# generate snp list yaml

mkdir -p full_traits_snp_list

wkdir=`pwd`
cd ../

for i in `seq 1 22`
do
  python generate_snp_list.py --yaml-template <(cat $wkdir/template.yaml | sed "s/XXX/$i/g") --trait-yaml $wkdir/full_traits.yaml --out-yaml $wkdir/full_traits_snp_list/chr$i.yaml
done
