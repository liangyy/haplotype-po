nsnps='1 10 50 100 500 1000'

if [[ -f run_test_speed.log ]]
then
  rm run_test_speed.log
fi

for n in $nsnps
do
  python test_speed.py --nsnp $n >> run_test_speed.log  
done	

