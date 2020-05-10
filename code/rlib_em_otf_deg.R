# input: yf1, ..., yfp and ym1, ..., ymp
# input: h1 and h2
# output: Pr(Z = 1)
# otf stands for on-the-fly
# NOT sex-specific
em_algorithm_otf_deg = function(yf, ym, h1, h2, tol = 1e-5, maxiter = 100, prior_z = 0.5, add_intercept = TRUE) {
  
  # add intercept
  if(add_intercept == TRUE) {
    h1 = cbind(rep(1, ncol(h1)), h1)
    h2 = cbind(rep(1, ncol(h2)), h2)
  }
  
  
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
  beta = matrix(0, ncol = p, nrow = k)
  sigma2 = matrix(1, ncol = p, nrow = 1)
  
  diff = tol + 1
  niter = 1
  lld = c()
  while(diff > tol & niter <= maxiter) {
    
    # E step
    l0 = calc_l(yf, ym, h1, h2, beta, sigma2, z = 0)
    l1 = calc_l(yf, ym, h1, h2, beta, sigma2, z = 1)
    gamma = calc_gamma(l0, l1, prior_z = prior_z)
    lld = c(lld, get_lld(l0, l1, sigma2, n, prior_z = prior_z))
    
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
    beta = update_beta(N, M, yf, tilde_N, tilde_M, ym)
    sigma2 = update_sigma2(n, yf, N, ym, tilde_N, beta)
    
    diff_b = calc_diff(beta_old, beta)
    diff_s = calc_diff(sigma2_old, sigma2)
    
    diff = diff_b + diff_s
    niter = niter + 1
    
  }
  
  # last update
  l0 = calc_l(yf, ym, h1, h2, beta, sigma2, z = 0)
  l1 = calc_l(yf, ym, h1, h2, beta, sigma2, z = 1)
  lld = c(lld, get_lld(l0, l1, sigma2, n, prior_z = prior_z))
  gamma = calc_gamma(l0, l1, prior_z = prior_z)
  
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
    yf_ = h1 %*% beta
    ym_ = h2 %*% beta
    s2f_ = sigma2
    s2m_ = sigma2
  } else if(z == 0) {
    ym_ = h1 %*% beta
    yf_ = h2 %*% beta
    s2m_ = sigma2
    s2f_ = sigma2
  }
  
  rf = yf - yf_  # n x p
  rm = ym - ym_
  
  l = - rowSums(mat_vec_div_by_row(rf ^ 2 / 2, s2f_)) - rowSums(mat_vec_div_by_row(rm ^ 2 / 2, s2m_))
  return(l)
  
}
calc_gamma = function(l0, l1, prior_z = 0.5) {
  max_l = pmax(l0, l1)
  r1 = l1 - max_l
  r0 = l0 - max_l
  o = exp(r1) * prior_z / (exp(r0) * (1 - prior_z) + exp(r1) * prior_z)
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
update_beta = function(N, M, y, N_, M_, y_) {
  beta = solve(N + N_, t(M) %*% y + t(M_) %*% y_)
  return(beta)
}
update_sigma2 = function(n, y, N, y_, N_, beta) {
  sigma2 = 1 / 2 / n * (diag(t(y) %*% y) - diag(t(beta) %*% N %*% beta) + diag(t(y_) %*% y_) - diag(t(beta) %*% N_ %*% beta))
  return(sigma2)
}
calc_diff = function(o, n) {
  out = sum((o - n) ^ 2)
  return(out)
}
logsum_ = function(a, b, prior_a = 0.5) {
  m = pmax(a, b)
  ma = a - m
  mb = b - m
  o = log(exp(ma) * prior_a + exp(mb) * (1 - prior_a)) + m
  return(o)
}
get_lld = function(l0, l1, sigma2, n, prior_z = 0.5) {
  l1_plus_l0 = logsum_(l0, l1, 1 - prior_z)
  o = sum(l1_plus_l0) - n / 2 * sum(log(sigma2) + log(sigma2)) - n * log(2 * pi)
  return(o)
}
