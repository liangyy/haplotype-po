# ARGS1: gene list
# ARGS2: chromosome number

genelist=$1

if [[ ! -d logs ]]
then
  mkdir -p logs
fi

for i in `cat $genelist`
do
  bash run_ld_clump_one_chromosome.screen $i $2
done

