em_per_snp = function(yf, ym, h1, h2, prior_prob_z = 0.5, maxiter = 1000, tol = 1e-10) {
  
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
  
  beta_f = matrix(0, nrow = k, ncol = p)
  inter_f = matrix(0, nrow = k, ncol = p)
  sigma2_f = matrix(1, nrow = k, ncol = p)
  beta_m = matrix(0, nrow = k, ncol = p)
  inter_m = matrix(0, nrow = k, ncol = p)
  sigma2_m = matrix(1, nrow = k, ncol = p)
  
  X_f = rbind(h1, h2)
  X_m = rbind(h2, h1)
  Y_f = rbind(yf, yf)
  Y_m = rbind(ym, ym)
  
  diff = tol + 1
  niter = 0
  lld = c()
  
  while(diff > tol & niter < maxiter) {
    
    # E step
    lf = get_l(yf, h1, h2, beta_f, inter_f, sigma2_f)
    lm = get_l(ym, h2, h1, beta_m, inter_m, sigma2_m)
    lld = c(lld, get_lld(lf, lm, prior_prob_z, sigma2_f, sigma2_m))
    gamma = get_gamma(lf, lm, prior_prob_z)
    
    # M step
    beta_f_old = beta_f
    beta_m_old = beta_m
    inter_f_old = inter_f
    inter_m_old = inter_m
    sigma2_f_old = sigma2_f
    sigma2_m_old = sigma2_m
    tmp_ = update_beta_and_inter(X_f, Y_f, gamma)
    beta_f = tmp_$beta
    inter_f = tmp_$inter
    sigma2_f = update_sigma2(X_f, Y_f, gamma, beta_f, inter_f)
    tmp_ = update_beta_and_inter(X_m, Y_m, gamma)
    beta_m = tmp_$beta
    inter_m = tmp_$inter
    sigma2_m = update_sigma2(X_m, Y_m, gamma, beta_m, inter_m)
    
    # diff
    diff_b = get_diff(beta_f_old, beta_f) + get_diff(beta_m_old, beta_m) 
    diff_i = get_diff(inter_f_old, inter_f) + get_diff(inter_m_old, inter_m) 
    diff_s = get_diff(sigma2_f_old, sigma2_f) + get_diff(sigma2_m_old, sigma2_m)
    diff = diff_s + diff_b + diff_i
    
    niter = niter + 1
  }
  # final update
  lf = get_l(yf, h1, h2, beta_f, inter_f, sigma2_f)
  lm = get_l(ym, h2, h1, beta_m, inter_m, sigma2_m)
  lld = c(lld, get_lld(lf, lm, prior_prob_z, sigma2_f, sigma2_m))
  gamma = get_gamma(lf, lm, prior_prob_z)
  
  list(prob_z = gamma, beta = list(f = beta_f, m = beta_m), inter = list(f = inter_f, m = inter_m), sigma2 = list(f = sigma2_f, m = sigma2_m), lld = lld)
}

get_lld = function(lf, lm, prior, sigma2f, sigma2m) {
  l1 = rowSums(lf$l1 + lm$l1) + log(prior)
  l0 = rowSums(lf$l0 + lm$l0) + log(1 - prior)
  lld = .logsum(l1, l0)
  o = sum(lld)
  o - nrow(lf$l1) / 2 * (sum(log(sigma2f)) + sum(log(sigma2m)))
}

get_l = function(y, hh1, hh2, beta, inter, sigma2) {
  l0 = matrix(0, ncol = ncol(y), nrow = nrow(y))
  l1 = matrix(0, ncol = ncol(y), nrow = nrow(y))
  # loop over trait
  for(pp in 1 : ncol(y)) {
    yp1 = .mat_vec_by_row(.mat_vec_by_row(hh1, beta[, pp], '*'), inter[, pp], '+')
    yp2 = .mat_vec_by_row(.mat_vec_by_row(hh2, beta[, pp], '*'), inter[, pp], '+')
    residual1 = .mat_vec_by_row(.mat_vec_by_col(yp1, y[, pp], '-') ^ 2, 2 * sigma2[, pp], '/') 
    residual2 = .mat_vec_by_row(.mat_vec_by_col(yp2, y[, pp], '-') ^ 2, 2 * sigma2[, pp], '/') 
    l0[, pp] = rowSums(residual2)
    l1[, pp] = rowSums(residual1)
  }
  list(l1 = -l1, l0 = -l0)
}

get_gamma = function(lf, lm, prior) {
  l1 = rowSums(lf$l1 + lm$l1) + log(prior)
  l0 = rowSums(lf$l0 + lm$l0) + log(1 - prior)
  lld = .logsum(l1, l0)
  exp(l1 - lld)
}

update_beta_and_inter = function(X, Y, gamma) {
  
  Gamma = c(gamma, 1 - gamma)
  # tilde_X = .mat_vec_by_col(X, sqrt(tilde_gamma), '*')
  # tilde_Y = .mat_vec_by_col(Y, sqrt(tilde_gamma), '*')
  
  w_bar = mean(Gamma)
  wx_bar = colMeans(.mat_vec_by_col(X, Gamma, '*'))
  wxsq_bar = colMeans(.mat_vec_by_col(X ^ 2, Gamma, '*'))
  
  # loop over traits
  inter = matrix(0, ncol = ncol(Y), nrow = ncol(X))
  beta = matrix(0, ncol = ncol(Y), nrow = ncol(X))
  
  for(pp in 1 : ncol(Y)) {
    wy_bar = mean(Y[, pp] * Gamma)
    wxy_bar = colMeans(.mat_vec_by_col(X, Y[, pp] * Gamma, '*'))
    inter[, pp] = ( (wxsq_bar * wy_bar) - (wx_bar * wxy_bar) ) / (wxsq_bar * w_bar - wx_bar ^ 2)
    beta[, pp] = ( - (wx_bar * wy_bar) + (w_bar * wxy_bar) ) / (wxsq_bar * w_bar - wx_bar ^ 2)
  }
  list(beta = beta, inter = inter)
}

update_sigma2 = function(X, Y, gamma, beta, inter) {
  
  tilde_gamma = c(gamma, 1 - gamma)
  tilde_ones = .mat_vec_by_col(matrix(1, ncol = ncol(X), nrow = nrow(X)), sqrt(tilde_gamma), '*')
  tilde_X = .mat_vec_by_col(X, sqrt(tilde_gamma), '*')
  tilde_Y = .mat_vec_by_col(Y, sqrt(tilde_gamma), '*')
  
  # loop over traits
  sigma2 = matrix(0, ncol = ncol(Y), nrow = ncol(X))
  for(pp in 1 : ncol(Y)) {
    Yp = .mat_vec_by_row(tilde_X, beta[, pp], '*') + .mat_vec_by_row(tilde_ones, inter[, pp], '*')
    R = .mat_vec_by_col(Yp, tilde_Y[, pp], '-') ^ 2
    sigma2[, pp] = colSums(R) / sum(tilde_gamma)
  }
  sigma2
}

get_diff = function(a, b) {
  sum((a - b) ^ 2)
}

.mat_vec_by_row = function(mat, vec, func) {
  sweep(mat, 2, vec, FUN = func)
}

.mat_vec_by_col = function(mat, vec, func) {
  sweep(mat, 1, vec, FUN = func)
}

.logsum = function(a, b) {
  # return log(exp(a) + exp(b))
  m = pmax(a, b)
  ma = a - m
  mb = b - m
  log(exp(ma) + exp(mb)) + m
}

.check_dim = function(mat, m, n) {
  if(ncol(mat) != n | nrow(mat) != m) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}