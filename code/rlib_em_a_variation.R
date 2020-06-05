# algorithm 1 (a variation)
# model:
#   y = b0 * covariates + b * g + epsilon  # this model is sex specific
#   with b >= 0
# note:
#   by "sex-specific", we mean b, b0, and the variance of epsilon are sex-specific
# 
# see note at ../analysis/init_idea_a_variation.R

em_algorithm_a_variation = function(y_father, y_mother, g_1, g_2, covar = NULL, maxiter = 15, tol = 1e-10, random_init = FALSE, non_negative = TRUE) {
  n = nrow(y_father)  # sample size
  p = ncol(y_father)  # number of phenotypes
  
  cc = matrix(1, ncol = 1, nrow = n)  # add intercept
  if(!is.null(covar)) {
    cc = cbind(cc, covar)  # add other covariates
  }
  
  Nc = ncol(cc)  # number of covariates
  
  sigma2_f = rep(1, p)
  sigma2_m = rep(1, p)
  # beta_f = rep(0, p)
  # beta_m = rep(0, p)
  beta_f = matrix(0, ncol = p, nrow = 1 + Nc)  # NOTE: PRS goes first!
  beta_m = matrix(0, ncol = p, nrow = 1 + Nc)
  if(random_init == TRUE) {
    beta_f = matrix(rnorm(p * (1 + Nc)), ncol = p, nrow = 1 + Nc)
    beta_m = matrix(rnorm(p * (1 + Nc)), ncol = p, nrow = 1 + Nc)
  }
  
  lld = c()
  niter = 0
  diff = tol + 1
  while(diff > tol & niter < maxiter) {
    
    # E step
    l1 = nlog_prob_y_given_z(
      y_father, y_mother, 
      g_1, g_2, cc, 
      beta_f, beta_m, sigma2_f, sigma2_m
    )
    l0 = nlog_prob_y_given_z(
      y_father, y_mother, 
      g_2, g_1, cc, 
      beta_f, beta_m, sigma2_f, sigma2_m
    )
    lld = c(lld, get_lld(l1, l0))  # lld = c(lld, -get_lld(l1, l0, sigma2_f, sigma2_m))
    gamma = 1 / (1 + exp(l1 - l0))
    
    # M step
    beta_f_old = beta_f; beta_m_old = beta_m; sigma2_f_old = sigma2_f; sigma2_m_old = sigma2_m
    beta_f = update_beta(y_father, g_1, g_2, cc, omega = c(gamma, 1 - gamma), non_negative)
    beta_m = update_beta(y_mother, g_2, g_1, cc, omega = c(gamma, 1 - gamma), non_negative)
    sigma2_f = update_sigma(y_father, g_1, g_2, beta_f, cc, omega = c(gamma, 1 - gamma))
    sigma2_m = update_sigma(y_mother, g_2, g_1, beta_m, cc, omega = c(gamma, 1 - gamma))
    
    # calc_diff
    diff_beta = get_diff(beta_f_old, beta_f) + get_diff(beta_m_old, beta_m) 
    diff_sigma2 = get_diff(sigma2_f_old, sigma2_f) + get_diff(sigma2_m_old, sigma2_m)
    diff = diff_beta + diff_sigma2
    
    # add on iteration
    niter = niter + 1
  }
  return(list(z_prob_n = gamma, sigma2_f = sigma2_f, sigma2_m = sigma2_m, beta_f = beta_f, beta_m = beta_m, lld = lld))
}
# END of MAIN algorithm 1

.mat_vec_by = function(mat, vec, func, margin, tag) {
  if(length(vec) != dim(mat)[margin]) {
    message('Wrong dimension in ', tag)
  }
  sweep(mat, margin, vec, FUN = func)
}

.mat_vec_by_row = function(mat, vec, func) {
  .mat_vec_by(mat, vec, func, 2, '.mat_vec_by_row')
  # if(length(vec) != ncol(mat)) {
  #   message('Wrong dimension in .mat_vec_by_row')
  # }
  # sweep(mat, 2, vec, FUN = func)
}

.mat_vec_by_col = function(mat, vec, func) {
  .mat_vec_by(mat, vec, func, 1, '.mat_vec_by_col')
  # if(length(vec) != nrow(mat)) {
  #   message('Wrong dimension in .mat_vec_by_col')
  # }
  # sweep(mat, 1, vec, FUN = func)
}

get_residual = function(y, g, covar, beta) {
  yg = .mat_vec_by_row(g, beta[1, ], '*')
  ycovar = covar %*% beta[-1, ]
  y - yg - ycovar
}

nlog_prob_y_given_z = function(
  yf, ym, 
  gf, gm, covar, 
  beta_f, beta_m, sigma2_f, sigma2_m
) {
  Nsample = nrow(yf)
  res_f = get_residual(yf, gf, covar, beta_f)
  res_m = get_residual(ym, gm, covar, beta_m)
  ratio_f = .mat_vec_by_row(res_f ^ 2, 2 * sigma2_f, '/')
  ratio_m = .mat_vec_by_row(res_m ^ 2, 2 * sigma2_m, '/')
  lf_n_by_p = .mat_vec_by_row(ratio_f, 1 / 2 * log(sigma2_f), '+')
  lm_n_by_p = .mat_vec_by_row(ratio_m, 1 / 2 * log(sigma2_m), '+')
  rowSums(lf_n_by_p) + rowSums(lm_n_by_p)
  # rowSums(ratio_f) + rowSums(ratio_m)
}

.logsum = function(lx, ly) {
  lmax = pmax(lx, ly)
  lx_ = lx - lmax
  ly_ = ly - lmax
  lmax + log(exp(lx_) + exp(ly_))
}

get_lld = function(l1, l0) {
  # sum(l1, l0)
  # sum(.logsum(-l1, -l0)) - length(l1) / 2 * sum(log(s2f) + log(s2m))
  sum(.logsum(-l1, -l0))
}

update_beta = function(y, g_1, g_2, cc, omega, non_negative) {
  
  # update all phenotypes simultaneously
  # Equation:
  #   Y: n x p
  #   X: n x p
  #   C: n x Nc
  #   A = (CtC)^-1 CtY  # Nc x p
  #   B = (CtC)^-1 CtX  # Nc x p
  #   D = CtY  # Nc x p
  #   E = XtY  # p x 1 (take diag)
  #   S = XtX - XtC B  # p x 1 (take diag)  
  #   BtD = Bt D  # p x 1 (take diag)
  #   beta = - S^-1 (BtD - E)  # p x 1
  #   beta_c = A - B beta (by row)  # Nc x p
  
  CtC = t(cc) %*% cc
  CtY = t(cc) %*% y
  Y = .mat_vec_by_col(rbind(y, y), sqrt(omega), '*')
  X = .mat_vec_by_col(rbind(g_1, g_2), sqrt(omega), '*')
  C = .mat_vec_by_col(rbind(cc, cc), sqrt(omega), '*')
  CtX = t(C) %*% X
  A = solve(CtC, CtY)
  B = solve(CtC, CtX)
  D = CtY
  E = colSums(X * Y)
  S = colSums(X * X) - colSums(CtX * B)
  BtD = colSums(B * D)
  beta = - (BtD - E) / S
  # the extra constraint on beta so that beta is non negative
  if(non_negative == TRUE) {
    beta[beta < 0] = 0
  }
  beta_c = A - .mat_vec_by_row(B, beta, '*')
  rbind(beta, beta_c)
}

update_sigma = function(y, g_1, g_2, beta, cc, omega) {
  r1 = get_residual(y, g_1, cc, beta)
  r2 = get_residual(y, g_2, cc, beta)
  r_tilde = .mat_vec_by_col(rbind(r1, r2), sqrt(omega), '*')
  denom = sum(omega)
  sigma2 = colSums(r_tilde ^ 2) / denom
  sigma2
}

get_diff = function(old, new) {
  sum((old - new) ^ 2)
}
