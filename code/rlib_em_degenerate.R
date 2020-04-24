# algorithm 0
# model:
#   y = intercept + g + epsilon  # this model is NOT sex-specific
# note:
#   by "NOT sex-specific", we mean the two sexes share the same intercept and variance 
#   intercept is called offset in the implementation
log_prob_y_g_given_z_deg = function(yf, ym, g1, g2, sigma2_p, offset, z) {
  # message('dim = ', dim(diff_y_g(yf, ym, g1, g2, z)))
  residuals = diff_y_g_deg(yf, ym, g1, g2, offset, z)
  o_father = - sweep(residuals$father, 2, FUN = '/',  2 * sigma2_p)
  o_mother = - sweep(residuals$mother, 2, FUN = '/',  2 * sigma2_p)
  o = o_father + o_mother
  # message('length = ', length(sigma2_p))
  return(rowSums(o))
}

subtract_offset_deg = function(xmat, offset) {
  return(
    sweep(xmat, 2, FUN = '-', offset)
  )
}

diff_y_g_deg = function(yf, ym, g1, g2, offset, z) {
  y_minus_g_ = y_minus_g_deg(yf, ym, g1, g2, z)
  o = list(
    father = subtract_offset_deg(y_minus_g_$father, offset) ^ 2,
    mother = subtract_offset_deg(y_minus_g_$mother, offset) ^ 2
  )
  return(o)
}

y_minus_g_deg = function(yf, ym, g1, g2, z) {
  if(z == 1) {
    o = list(father = yf - g1, mother = ym - g2)
  } else if(z == 0) {
    o = list(mother = ym - g1, father = yf - g2)
  } 
  return(o)
}

em_algorithm_deg = function(y_father, y_mother, g_1, g_2, maxiter = 15, prior_z = 0.5) {
  n = nrow(y_father)  # sample size
  p = ncol(y_father)  # number of phenotypes
  # z_prob_n = rep(0.5, n)
  sigma2_p = rep(1, p)
  offset = rep(0, p)
  lld = c()
  niter = 0
  while(niter < maxiter) {
    
    # E step
    ## to handle potential over/underflow
    lp1 = log_prob_y_g_given_z_deg(
      y_father, y_mother, 
      g_1, g_2, 
      sigma2_p, 
      offset, 
      1
    )
    lp0 = log_prob_y_g_given_z_deg(
      y_father, y_mother, 
      g_1, g_2, 
      sigma2_p, 
      offset, 
      0
    )
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
    y_minus_g_1 = y_minus_g_deg(y_father, y_mother, g_1, g_2, 1)
    y_minus_g_0 = y_minus_g_deg(y_father, y_mother, g_1, g_2, 0)
    offset1 = 1 / n * colSums(y_minus_g_1$father * w + y_minus_g_0$father * (1 - w))
    offset2 = 1 / n * colSums(y_minus_g_1$mother * w + y_minus_g_0$mother * (1 - w))
    offset = (offset1 + offset2) / 2
    residual_1 = diff_y_g_deg(y_father, y_mother, g_1, g_2, offset, 1)
    residual_0 = diff_y_g_deg(y_father, y_mother, g_1, g_2, offset, 0)
    sigma2_p1 = 1 / (n - 1) * colSums(sweep(residual_1$father, 1, FUN = '*', w)  + sweep(residual_0$father, 1, FUN = '*', 1 - w))
    sigma2_p2 = 1 / (n - 1) * colSums(sweep(residual_1$mother, 1, FUN = '*', w)  + sweep(residual_0$mother, 1, FUN = '*', 1 - w))
    sigma2_p = (sigma2_p1 + sigma2_p2) / 2
    
    
    # some others
    niter = niter + 1
  }
  return(list(z_prob_n = w, sigma2_p = sigma2_p, offset = offset, lld = lld))
}
# END of algorithm 0
