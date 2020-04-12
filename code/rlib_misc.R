unfold_parent_illness_code = function(dd, tbl_code, tbl_name) {
  out = matrix(0, nrow = nrow(dd), ncol = length(tbl_code))
  for(i in 1 : ncol(dd)) {
    dcol = dd[, i]
    for(j in 1 : ncol(out)) {
      has_code = !is.na(dcol) & dcol == tbl_code[j]
      out[has_code, j] = out[has_code, j] + 1
    }
  }
  out = as.data.frame(out)
  colnames(out) = tbl_name
  return(out)
}

# use to clean up the table copy/pasted from data field Data tab
clean_up = function(filename) {
  e = read.delim2(filename)
  write.table(e, filename, sep = ',', col.names = T, row.names = F, quote = F)
}