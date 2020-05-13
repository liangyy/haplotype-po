# implement logistic EM
# Inputs:
#   yf: father phenotype
#   ym: mother phenotype
#   h1: haplotype 1
#   h2: haplotype 2
#   prior prob z: prior_prob_z
# Outputs:
#   prob_z

em_logistic = function(yf, ym, h1, h2, prior_prob_z, tol = 1e-10, maxiter = 100) {
  
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
  
  X_f = rbind(h1, h2)
  X_m = rbind(h2, h1)
  Y_f = rbind(yf, yf)
  Y_m = rbind(ym, ym)
  beta_f = matrix(0, ncol = p, nrow = k)
  beta_m = matrix(0, ncol = p, nrow = k)
  lld = c()
  
  diff = tol + 1
  niter = 0
  while(diff > tol & niter < maxiter) {
    mu_f = get_mu(beta_f, X_f)
    mu_m = get_mu(beta_m, X_m)
    
    # E step
    l_f = get_l(mu_f, Y_f)
    l_m = get_l(mu_m, Y_m)
    gamma = get_gamma(l_f, l_m, prior_prob_z)
    lld = c(lld, sum(get_lld(l_f, l_m, prior_prob_z)))
    
    # archive
    beta_f_old = beta_f
    beta_m_old = beta_m
    
    # M step
    beta_f = solve_weighted_logistic(X_f, Y_f, c(gamma, 1 - gamma), beta_init = beta_f)
    beta_m = solve_weighted_logistic(X_m, Y_m, c(gamma, 1 - gamma), beta_init = beta_m)
    
    # diff
    diff_f = get_diff(beta_f_old, beta_f)
    diff_m = get_diff(beta_m_old, beta_m)
    diff = diff_f + diff_m
    
    niter = niter + 1
  }
  
  # final calculation
  mu_f = get_mu(beta_f, X_f)
  mu_m = get_mu(beta_m, X_m)
  l_f = get_l(mu_f, Y_f)
  l_m = get_l(mu_m, Y_m)
  gamma = get_gamma(l_f, l_m, prior_prob_z)
  lld = c(lld, sum(get_lld(l_f, l_m, prior_prob_z)))
  return(list(beta_f = beta_f, beta_m = beta_m, prob_z = gamma, lld = lld))
}

solve_weighted_logistic = function(X, Y, gamma, beta_init) {
  p = dim(Y)[2]
  k = dim(beta_init)[1]
  beta = matrix(0, nrow = k, ncol = p)
  for(pp in 1 : p) {
    beta[, pp] = .solve_weighted_logstic(X, Y[, pp], gamma, beta_init[, pp])
  }
  beta
}

get_diff = function(a, b) {
  mean((a - b) ^ 2)
}

get_lld = function(lf, lm, prior) {
  l1 = rowSums(log(lf$l1) + log(lm$l1) + log(prior))
  l0 = rowSums(log(lf$l0) + log(lm$l0) + log(1 - prior))
  .logsum(l1, l0)
}

get_gamma = function(lf, lm, prior) {
  l1 = rowSums(log(lf$l1) + log(lm$l1) + log(prior))
  lld = get_lld(lf, lm, prior)
  exp(l1) / exp(lld)
}

get_l = function(mu, y) {
  o = mu
  o[y == 0] = 1 - mu[y == 0]
  n = dim(o)[1] / 2
  list(l1 = o[1 : n, ], l0 = o[(n+1):(2*n), ])
}

get_mu = function(beta, X) {
  .logistic(X %*% beta)
}

.solve_weighted_logstic = function(X, y, gamma, beta_init, tol = 1e-10, maxiter = 100) {
  beta = beta_init
  
  sqrt_gamma = sqrt(gamma)
  tilde_X = .mat_vec_multiply_by_col(X, sqrt_gamma)
  tilde_y = sqrt_gamma * y
  Xy = t(tilde_X) %*% tilde_y
  
  diff = tol + 1
  niter = 0
  while(diff > tol & niter < maxiter) {
    mu = get_mu(beta, X)
    S = mu * (1 - mu)
    tilde_mu = sqrt_gamma * mu
    SX = .mat_vec_multiply_by_col(tilde_X, S)
    XSX = t(tilde_X) %*% SX
    XSX_beta = XSX %*% beta
    XMu = t(tilde_X) %*% tilde_mu
    
    # archive
    beta_old = beta
    
    # lin solve
    LHS = XSX
    RHS = XSX_beta + Xy - XMu
    beta = solve(LHS, RHS)
    
    # diff
    diff = get_diff(beta, beta_old)
    
    niter = niter + 1
  }
  beta
}

.logsum = function(a, b) {
  # return log(exp(a) + exp(b))
  m = pmax(a, b)
  ma = a - m
  mb = b - m
  log(exp(ma) + exp(mb)) + m
}

.logistic = function(u) {
  1 / (1 + exp(-u))
}

.mat_vec_multiply_by_col = function(mat, vec) {
  sweep(mat, 1, vec, FUN = '*')
}

.check_dim = function(mat, m, n) {
  if(ncol(mat) != n | nrow(mat) != m) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}