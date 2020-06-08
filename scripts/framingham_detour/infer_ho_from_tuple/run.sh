# ARGS1: config file 
# ARGS2: log file
# ARGS3: number of cores 

$SNMK --configfile $1 -s infer.snmk --cores $3 -p > $2 2>&1

