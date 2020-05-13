n = 100
k = 5
p = 6
h1 = matrix(rnorm(n * k), ncol = k)
h2 = matrix(rnorm(n * k), ncol = k)
beta_f = matrix(rnorm(k * p), ncol = p)
beta_m = matrix(rnorm(k * p), ncol = p)
y_f = matrix(0, ncol = p, nrow = n)
y_m = matrix(0, ncol = p, nrow = n)
u_f = h1 %*% beta_f
u_m = h2 %*% beta_m

for(pp in 1 : p) {
  y_f[, pp] = rbinom(n, size = 1, prob = .logistic(u_f[, pp]))
  y_m[, pp] = rbinom(n, size = 1, prob = .logistic(u_m[, pp]))
}

o = em_logistic(y_f, y_m, h1, h2, prior_prob_z = 0.5)
hist(o$prob_z)
