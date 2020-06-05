y = matrix(rnorm(100), ncol = 2)
g1 = matrix(rnorm(100), ncol = 2)
g2 = matrix(rnorm(100), ncol = 2)
co = matrix(rnorm(150), ncol = 3)
gamma_ = runif(50)
beta = update_beta(y, g1, g2, co, omega = c(gamma_, 1 - gamma_))
update_sigma(y, g1, g2, beta, co, omega = c(gamma_, 1 - gamma_))


test_solver = function(y, g1, g2, co, omega) {
  y = c(y, y) 
  g = c(g1, g2) 
  co = rbind(co, co)
  x = cbind(g, co)
  lm(y ~ -1 + x, weights = omega)
}

test_solver(y[, 1], g1[, 1], g2[, 1], co, omega = c(gamma_, 1 - gamma_))
test_solver(y[, 2], g1[, 2], g2[, 2], co, omega = c(gamma_, 1 - gamma_))


test_sigma2 = function(y, g1, g2, co, beta, beta0, omega) {
  r1 = y - g1 * beta - co %*% beta0
  r2 = y - g2 * beta - co %*% beta0
  sum(c(r1 ^ 2, r2 ^ 2) * omega) / sum(omega)
}

test_sigma2(y[, 1], g1[, 1], g2[, 1], co, beta[1, 1], beta[-1, 1], omega = c(gamma_, 1 - gamma_))
test_sigma2(y[, 2], g1[, 2], g2[, 2], co, beta[1, 2], beta[-1, 2], omega = c(gamma_, 1 - gamma_))
