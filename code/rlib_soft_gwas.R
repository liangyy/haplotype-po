soft_gwas = function(y, x1, x2, weights, niter = 10000) {
  n = ncol(x1)
  intercept = matrix(1, nrow = nrow(x1), ncol = 1)
  blist = c()
  selist = c()
  lrtlist = c()
  l1list = c()
  l0list = c()
  s2list = c()
  for(i in 1 : n) {
    xx1 = cbind(intercept, x1[, i])
    xx2 = cbind(intercept, x2[, i])
    yy = y[, i, drop = FALSE]
    ww = weights[, i]
    o = em_algorithm_for_gwas(yy, xx1, xx2, prior_z = ww, maxiter = niter, tol = 1e-10)
    o0 = em_algorithm_for_gwas(yy, intercept, intercept, prior_z = ww, maxiter = niter, tol = 1e-10)
    l1 = o$lld[length(o$lld)]
    l0 = o0$lld[length(o0$lld)]
    lrt_stat = 2 * (l1 - l0)
    
    ## FIXME
    if(lrt_stat < 0) {
      lrt_stat = 0
    }
    ## END
    
    pval = exp(pchisq(lrt_stat, 1, lower.tail = FALSE, log.p = TRUE))
    zval = sqrt(lrt_stat)
    bhat = o$beta[2] 
    bhat_se = abs(bhat / zval)
    blist = c(blist, bhat)
    selist = c(selist, bhat_se)
    lrtlist = c(lrtlist, lrt_stat)
    l1list = c(l1list, l1)
    l0list = c(l0list, l0)
    s2list = c(s2list, o$sigma2)
  }
  return(list(bhat = blist, bhat_se = selist, lrt_stat = lrtlist, l1 = l1list, l0 = l0list, sigma2 = s2list))
}

.calc_l = function(y, x, b, s2) {
  r = y - x %*% b
  ll = - r ^ 2 / 2 / s2 - 1 / 2 * log(s2)
  ll
}

.logsum = function(a, b, prior_a = 0.5) {
  m = pmax(a, b)
  ma = a - m
  mb = b - m
  o = log(exp(ma) * prior_a + exp(mb) * (1 - prior_a)) + m
  return(o)
}

.get_lld = function(l0, l1, prior_z = 0.5) {
  l1_plus_l0 = .logsum(l0, l1, 1 - prior_z)
  o = sum(l1_plus_l0) 
  return(o)
}

.calc_gamma = function(l0, l1, prior_z = 0.5) {
  max_l = pmax(l0, l1)
  r1 = l1 - max_l
  r0 = l0 - max_l
  o = exp(r1) * prior_z / (exp(r0) * (1 - prior_z) + exp(r1) * prior_z)
  return(o)
}

.update_beta = function(y, x1, x2, w) {
  ww = sqrt(c(w, 1 - w))
  yy = c(y, y) * ww
  xx = sweep(rbind(x1, x2), 1, ww, FUN = '*')
  xtx = t(xx) %*% xx
  xty = t(xx) %*% yy
  beta = solve(xtx, xty)
  beta
}

.update_sigma2 = function(y, x1, x2, gamma, beta) {
  r1 = y - x1 %*% beta
  r2 = y - x2 %*% beta
  sigma2 = mean(gamma * r1 ^ 2 + (1 - gamma) * r2 ^ 2)
  sigma2
}

.calc_diff = function(new, old) {
  sum((new - old) ^ 2)
}

em_algorithm_for_gwas = function(y, x1, x2, prior_z, maxiter = 1000, tol = 1e-10) {
  beta = rep(0, ncol(x1))
  sigma2 = 1
  
  diff = tol + 1
  niter = 1
  lld = c()
  while(diff > tol & niter <= maxiter) {
    # E step
    l0 = .calc_l(y, x2, beta, sigma2)
    l1 = .calc_l(y, x1, beta, sigma2)
    lld = c(lld, .get_lld(l0, l1, prior_z = prior_z))
    gamma = .calc_gamma(l0, l1, prior_z = prior_z)
    # M step
    beta_old = beta
    sigma2_old = sigma2
    beta = .update_beta(y, x1, x2, gamma)
    sigma2 = .update_sigma2(y, x1, x2, gamma, beta)
    
    diff_b = .calc_diff(beta_old, beta)
    diff_s = .calc_diff(sigma2_old, sigma2)
    diff = diff_b + diff_s
    
    niter = niter + 1
  }
  # just one more E step
  l0 = .calc_l(y, x2, beta, sigma2)
  l1 = .calc_l(y, x1, beta, sigma2)
  lld = c(lld, .get_lld(l0, l1, prior_z = prior_z))
  gamma = .calc_gamma(l0, l1, prior_z = prior_z)
  
  return(list(beta = beta, sigma2 = sigma2, lld = lld, gamma = gamma))
}
