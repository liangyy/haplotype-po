# ARGS1: gene list
# ARGS2: chromosome number

genelist=$1

if [[ ! -d logs ]]
then
  mkdir -p logs
fi

for i in `cat $genelist`
do
  export TRAIT=$i
  export CHRNUM=$2
  bash run_ld_clump_one_chromosome.screen
done

