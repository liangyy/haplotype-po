# algorithm 1
# model:
#   y = g + epsilon  # this model shared by both parents
log_prob_y_g_given_z = function(yf, ym, g1, g2, sigma2_p, z) {
  # message('dim = ', dim(diff_y_g(yf, ym, g1, g2, z)))
  o = - sweep(diff_y_g(yf, ym, g1, g2, z), 2, FUN = '/',  2 * sigma2_p)
  # message('length = ', length(sigma2_p))
  return(rowSums(o))
}
diff_y_g = function(yf, ym, g1, g2, z) {
  if(z == 1) {
    o = (yf - g1) ^ 2 + (ym - g2) ^ 2
  } else if(z == 0) {
    o = (ym - g1) ^ 2 + (yf - g2) ^ 2
  } 
  return(o)
}
em_algorithm = function(y_father, y_mother, g_1, g_2, maxiter = 15, prior_z = 0.5) {
  n = nrow(y_father)  # sample size
  p = ncol(y_father)  # number of phenotypes
  # z_prob_n = rep(0.5, n)
  sigma2_p = rep(1, p)
  lld = c()
  niter = 0
  while(niter < maxiter) {
    
    # E step
    ## to handle potential over/underflow
    lp1 = log_prob_y_g_given_z(y_father, y_mother, g_1, g_2, sigma2_p, 1)
    lp0 = log_prob_y_g_given_z(y_father, y_mother, g_1, g_2, sigma2_p, 0)
    # message('lp1 = ', lp1[2], '; lp0 = ', lp0[2])
    lp_max = pmax(lp1, lp0)
    lp1 = lp1 - lp_max
    lp0 = lp0 - lp_max
    p1 = exp(lp1) * prior_z 
    p0 = exp(lp0) * (1 - prior_z)
    ## 
    w = p1 / (p1 + p0)
    lld = c(lld, sum(log((p1 + p0)) + lp_max))
    # message(p1[1], ';', p0[1], ' ;', w[1])
    
    # M step
    sigma2_p = 1 / 2 / n * colSums(sweep(diff_y_g(y_father, y_mother, g_1, g_2, 1), 1, FUN = '*', w)  + sweep(diff_y_g(y_father, y_mother, g_1, g_2, 0), 1, FUN = '*', 1 - w))
    # message('sigma2 = ', paste(sigma2_p, collapse = '; '))
    
    # some others
    niter = niter + 1
  }
  return(list(z_prob_n = w, sigma2_p = sigma2_p, lld = lld))
}
# END of algorithm 1
