# input: yf1, ..., yfp and ym1, ..., ymp
# input: h1 and h2
# output: Pr(Z = 1)
# otf stands for on-the-fly
em_algorithm_otf = function(yf, ym, h1, h2, tol = 1e-5, maxiter = 10) {
  
  # add intercept
  h1 = cbind(rep(1, ncol(h1)), h1)
  h2 = cbind(rep(1, ncol(h2)), h2)
  
  # dimensions
  p = ncol(yf)
  n = nrow(yf)
  k = ncol(h1)
  if(
    .check_dim(h1, n, k) &
    .check_dim(h2, n, k) &
    .check_dim(yf, n, p) &
    .check_dim(ym, n, p) == FALSE
  ) {
    message('Wrong dimension')
    return(NULL)
  }
  
  # initialize
  beta = list()
  sigma2 = list()
  beta$f = matrix(0, ncol = p, nrow = k)
  beta$m = matrix(0, ncol = p, nrow = k)
  sigma2$f = matrix(1, ncol = p, nrow = 1)
  sigma2$m = matrix(1, ncol = p, nrow = 1)
  
  diff = tol + 1
  niter = 1
  lld = c()
  while(diff > tol & niter <= maxiter) {
    
    # E step
    l0 = calc_l(yf, ym, h1, h2, beta, sigma2, z = 0)
    l1 = calc_l(yf, ym, h1, h2, beta, sigma2, z = 1)
    lld = c(lld, get_lld(l0, l1, sigma2, n))
    gamma = calc_gamma(l0, l1)
    
    # M step
    
    ## prepare M, N
    d_gamma_h1 = mat_vec_mul_by_col(h1, gamma)
    d_ngamma_h2 = mat_vec_mul_by_col(h2, 1 - gamma)
    M = d_gamma_h1 + d_ngamma_h2
    N = t(h1) %*% d_gamma_h1 + t(h2) %*% d_ngamma_h2
    tilde_M = h1 + h2 - M
    tilde_N = t(h1) %*% h1 + t(h2) %*% h2 - N
    
    ## update 
    beta_old = beta
    sigma2_old = sigma2
    beta$f = update_beta(N, M, yf)
    beta$m = update_beta(tilde_N, tilde_M, ym)
    sigma2$f = update_sigma2(n, yf, N, beta$f)
    sigma2$m = update_sigma2(n, ym, tilde_N, beta$m)
    
    diff_b = calc_diff(beta_old, beta)
    diff_s = calc_diff(sigma2_old, sigma2)
    
    diff = diff_b + diff_s
    niter = niter + 1
    
  }
  
  # last update
  l0 = calc_l(yf, ym, h1, h2, beta, sigma2, z = 0)
  l1 = calc_l(yf, ym, h1, h2, beta, sigma2, z = 1)
  lld = c(lld, get_lld(l0, l1, sigma2, n))
  gamma = calc_gamma(l0, l1)
  
  return(list(beta = beta, sigma2 = sigma2, gamma = gamma, lld = lld))
  
}
.check_dim = function(mat, m, n) {
  if(ncol(mat) != n | nrow(mat) != m) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
calc_l = function(yf, ym, h1, h2, beta, sigma2, z) {
  
  if(z == 1) {
    yf_ = h1 %*% beta$f
    ym_ = h2 %*% beta$m
    s2f_ = sigma2$f
    s2m_ = sigma2$m
  } else if(z == 0) {
    ym_ = h1 %*% beta$m
    yf_ = h2 %*% beta$f
    s2m_ = sigma2$m
    s2f_ = sigma2$f
  }
  
  rf = yf - yf_  # n x p
  rm = ym - ym_
  
  l = - rowSums(mat_vec_div_by_row(rf ^ 2 / 2, s2f_)) - rowSums(mat_vec_div_by_row(rm ^ 2 / 2, s2m_))
  return(l)
  
}
calc_gamma = function(l0, l1) {
  max_l = pmax(l0, l1)
  r1 = l1 - max_l
  r0 = l0 - max_l
  o = exp(r1) / (exp(r0) + exp(r1))
  return(o)
}
mat_vec_mul_by_col = function(mat, vec) {
  o = sweep(mat, 1, vec, FUN = '*')
  return(o)
}
mat_vec_div_by_row = function(mat, vec) {
  o = sweep(mat, 2, vec, FUN = '/')
  return(o)
}
update_beta = function(N, M, y) {
  beta = solve(N, t(M) %*% y)
  return(beta)
}
update_sigma2 = function(n, y, N, beta) {
  sigma2 = 1 / n * (diag(t(y) %*% y) - diag(t(beta) %*% N %*% beta))
  return(sigma2)
}
calc_diff = function(o, n) {
  out = 0
  for(i in names(o)) {
    out = out + sum((o[[i]] - n[[i]]) ^ 2)
  }
  return(out)
}
logsum_ = function(a, b) {
  m = pmax(a, b)
  ma = a - m
  mb = b - m
  o = log(exp(ma) + exp(mb)) + m
  return(o)
}
get_lld = function(l0, l1, sigma2, n) {
  l1_plus_l0 = logsum_(l0, l1)
  o = sum(l1_plus_l0) - n / 2 * sum(log(sigma2$f) + log(sigma2$m))
  return(o)
}