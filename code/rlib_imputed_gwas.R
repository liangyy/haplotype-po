imputed_gwas = function(y, x1, x2, weights) {
  n = ncol(x1)
  
  xx = x1 * weights + x2 * (1 - weights)
  
  x0 = matrix(1, nrow = nrow(xx), ncol = ncol(xx))
  x0 = sweep(x0, 2, FUN = '/', colSums(x0))
  
  a_y = apply(y, 2, sd)
  b_x = apply(xx, 2, sd)
  y_tilde = sweep(y, 2, FUN = '/', a_y) 
  x_tilde = sweep(xx, 2, FUN = '/', b_x) 
  xtx = colSums(x_tilde ^ 2)
  xbar = colSums(x_tilde * x0) 
  ybar = colSums(y_tilde * x0) 
  n = colSums(x0 ^ 2)
  xty = colSums(x_tilde * y_tilde)
  denom = n * xtx - xbar ^ 2
  mu0 = 1 / denom * (xtx * ybar - xbar * xty)
  bhat = 1 / denom * (- xbar * ybar + xty * n)
  ypred = sweep(x_tilde, 2, FUN = '*', bhat) + sweep(x0, 2, FUN = '*', mu0)
  sigma2 = 1 / (nrow(x_tilde) - 2) * colSums((ypred - y_tilde) ^ 2)
  bhat_se = sqrt(sigma2 * n / denom)
  # print(denom)
  bhat = bhat * a_y / b_x
  bhat_se = bhat_se * a_y / b_x
  return(list(bhat = bhat, bhat_se = bhat_se))
}