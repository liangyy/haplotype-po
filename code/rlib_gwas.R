# X: matrix n x p
# y: vector p x 1
# run y ~ Xj
run_gwas = function(X, y) {
  a_y = sd(y)
  b_x = apply(X, 2, sd)
  y_tilde = y / a_y
  # message(1)
  x_tilde = sweep(X, 2, FUN = '/', b_x)
  # message(2)
  xtx = colSums(x_tilde ^ 2)
  xbar = colMeans(x_tilde)
  ybar = mean(y_tilde)
  n = nrow(x_tilde)
  # message(dim(x_tilde), '  ', length(y_tilde))
  xty = colSums(sweep(x_tilde, 1, FUN = '*', y_tilde))
  # message(3)
  denom = xtx - n * xbar ^ 2
  mu0 = 1 / denom * (xtx * ybar - xbar * xty)
  bhat = 1 / denom * (-n * xbar * ybar + xty)
  ypred = sweep(sweep(x_tilde, 2, FUN = '*', bhat), 2, FUN = '+', mu0)
  # message(4)
  sigma2 = 1 / (n - 2) * colSums((sweep(ypred, 1, FUN = '-', y_tilde)) ^ 2)
  # message(5)
  bhat_se = sqrt(sigma2 / (xtx - n * xbar ^ 2))
  bhat = bhat * a_y / b_x
  bhat_se = bhat_se * a_y / b_x
  return(list(bhat = bhat, bhat_se = bhat_se))
}

# X: matrix n x p
# y: matrix n x p
# w: matrix n x p
# run yj ~ Xj with weights wj
run_gwas_pairwise = function(X, y, weights = NULL) {
  x0 = matrix(1, nrow = nrow(X), ncol = ncol(X))
  if(!is.null(weights)) {
    X = sqrt(weights) * X
    y = sqrt(weights) * y
    x0 = sqrt(weights) * x0
  }
  x0 = sweep(x0, 2, FUN = '/', colSums(x0))
  
  a_y = apply(y, 2, sd)
  b_x = apply(X, 2, sd)
  y_tilde = sweep(y, 2, FUN = '/', a_y) 
  x_tilde = sweep(X, 2, FUN = '/', b_x) 
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

tval2pval = function(tval, df) {
  exp(pt(abs(tval), df = df, lower.tail = FALSE, log.p = TRUE)) * 2
}

# # X: matrix n x p
# # y: matrix n x p
# # w: matrix n x p
# # run yj ~ Xj
# run_gwas_pairwise = function(X, y, weights = NULL) {
#   if(is.null(weights)) {
#     return(.run_gwas_pairwise(X, y))
#   } else {
#     X = sqrt(weights) * X
#     y = sqrt(weights) * y
#     return(.run_gwas_pairwise(X, y))
#   }
# }