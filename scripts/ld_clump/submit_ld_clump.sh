# ARGS1: gene list

genelist=$1

if [[ ! -d logs ]]
then
  mkdir -p logs
fi

for i in `cat $genelist`
do
  export TRAIT=$i
  bash run_ld_clump.screen
done

