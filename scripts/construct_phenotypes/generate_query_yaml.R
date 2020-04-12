e = "
# ARGS1: input query header including sample filters and covariates.
# ARGS2: a table of phenotype you'd like to query. 
#        Contain two columns: code (Data-Field ID), name (string ID as label).
#        In CSV format.
# ARGS2: output query yaml.
"

args = commandArgs(trailingOnly = TRUE)
if(args[1] == '-h') {
  message(e)
  quit()
}

library(dplyr)

current_yaml = yaml::read_yaml(args[1])

dict = read.csv("http://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv")

pheno = read.csv(args[2])

pheno = pheno %>% inner_join(dict %>% select(FieldID, Instances, Array), by = c('code' = 'FieldID'))

for(i in 1 : nrow(pheno)) {
  if(pheno$name[i] %in% names(current_yaml)) {
    message('Duplicated name: ', pheno$name[i], '. Exit!')
    quit()
  }
  pheno_name = as.character(pheno$name[i])
  pheno_code = pheno$code[i]
  current_yaml[[pheno_name]] = list()
  for(inst in 0 : (pheno$Instances[i] - 1)) {
    for(arr in 0 : (pheno$Array[i] - 1)) {
      name_ = paste0(pheno_name, '_', inst, '_', arr)
      id_ = paste0('c', pheno_code, '_', inst, '_', arr)
      current_yaml[[pheno_name]][[name_]] = id_
    }
  }
}  

yaml::write_yaml(current_yaml, args[3])
