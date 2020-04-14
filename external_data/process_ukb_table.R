source('code/rlib_misc.R')
for(i in c('father', 'mother')) {
  for(j in c(1, 2)) {
    clean_up(paste0('external_data/', i, '_illness_data_table_instance', j, '.csv'))
  }
}