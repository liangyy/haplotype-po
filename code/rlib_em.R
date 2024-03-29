# algorithm 1
# model:
#   y = intercept + g + epsilon  # this model is sex specific
# note:
#   by "sex-specific", we mean intercept, the variance of epsilon are both sex-specific
#   intercept is called offset in the implementation
log_prob_y_g_given_z = function(yf, ym, g1, g2, sigma2_p1, sigma2_p2, offset1, offset2, z) {
  # message('dim = ', dim(diff_y_g(yf, ym, g1, g2, z)))
  residuals = diff_y_g(yf, ym, g1, g2, offset1, offset2, z)
  o_father = - sweep(residuals$father, 2, FUN = '/',  2 * sigma2_p1)
  o_mother = - sweep(residuals$mother, 2, FUN = '/',  2 * sigma2_p2)
  o = o_father + o_mother
  # message('length = ', length(sigma2_p))
  return(rowSums(o))
}

subtract_offset = function(xmat, offset) {
  return(
    sweep(xmat, 2, FUN = '-', offset)
  )
}

diff_y_g = function(yf, ym, g1, g2, offset1, offset2, z) {
  y_minus_g_ = y_minus_g(yf, ym, g1, g2, z)
  o = list(
    father = subtract_offset(y_minus_g_$father, offset1) ^ 2,
    mother = subtract_offset(y_minus_g_$mother, offset2) ^ 2
  )
  return(o)
}

y_minus_g = function(yf, ym, g1, g2, z) {
  if(z == 1) {
    o = list(father = yf - g1, mother = ym - g2)
  } else if(z == 0) {
    o = list(mother = ym - g1, father = yf - g2)
  } 
  return(o)
}

em_algorithm = function(y_father, y_mother, g_1, g_2, maxiter = 15, prior_z = 0.5) {
  n = nrow(y_father)  # sample size
  p = ncol(y_father)  # number of phenotypes
  # z_prob_n = rep(0.5, n)
  sigma2_p1 = rep(1, p)
  sigma2_p2 = rep(1, p)
  offset1 = rep(0, p)
  offset2 = rep(0, p)
  lld = c()
  niter = 0
  while(niter < maxiter) {
    
    # E step
    ## to handle potential over/underflow
    lp1 = log_prob_y_g_given_z(
      y_father, y_mother, 
      g_1, g_2, 
      sigma2_p1, sigma2_p2, 
      offset1, offset2, 
      1
    )
    lp0 = log_prob_y_g_given_z(
      y_father, y_mother, 
      g_1, g_2, 
      sigma2_p1, sigma2_p2, 
      offset1, offset2, 
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
    y_minus_g_1 = y_minus_g(y_father, y_mother, g_1, g_2, 1)
    y_minus_g_0 = y_minus_g(y_father, y_mother, g_1, g_2, 0)
    offset1 = 1 / n * colSums(y_minus_g_1$father * w + y_minus_g_0$father * (1 - w))
    offset2 = 1 / n * colSums(y_minus_g_1$mother * w + y_minus_g_0$mother * (1 - w))
    residual_1 = diff_y_g(y_father, y_mother, g_1, g_2, offset1, offset2, 1)
    residual_0 = diff_y_g(y_father, y_mother, g_1, g_2, offset1, offset2, 0)
    sigma2_p1 = 1 / (n - 1) * colSums(sweep(residual_1$father, 1, FUN = '*', w)  + sweep(residual_0$father, 1, FUN = '*', 1 - w))
    sigma2_p2 = 1 / (n - 1) * colSums(sweep(residual_1$mother, 1, FUN = '*', w)  + sweep(residual_0$mother, 1, FUN = '*', 1 - w))

    
    # some others
    niter = niter + 1
  }
  return(list(z_prob_n = w, sigma2_p1 = sigma2_p1, sigma2_p2 = sigma2_p2, offset1 = offset1, offset2 = offset2, lld = lld))
}
# END of algorithm 1
