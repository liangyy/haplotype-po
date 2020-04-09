get_maf = function(p, maf_low, maf_high) {
  return(runif(p, min = maf_low, max = maf_high))
}
sim_hap = function(n, p, maf) {
  hlist = list()
  for(i in 1:2) {
    hlist[[i]] = matrix(rbinom(p * n, 1, rep(maf, n)), byrow = TRUE, nrow = n, ncol = p)
  }
  return(hlist)
}
transmit_haplo = function(father, mother) {
  out = list(
    matrix(NA, nrow = nrow(father[[1]]), ncol = ncol(father[[1]])),
    matrix(NA, nrow = nrow(father[[2]]), ncol = ncol(father[[2]]))
  )
  for(i in 1 : nrow(out[[1]])) {
    f_idx = sample(1:2, 1)
    m_idx = sample(1:2, 1)
    out[[1]][i, ] = father[[f_idx]][i, ]
    out[[2]][i, ] = mother[[m_idx]][i, ]
  }
  return(out)
}
spike_and_slab = function(n, pi_, sigma2) {
  b = rnorm(n, sd = sqrt(sigma2))
  b[runif(n) <= pi_] = 0
  return(b)
}
simulate_pheno = function(h, beta, h2, maf) {
  n = nrow(h[[1]])
  p = ncol(h[[1]])
  k = ncol(beta)
  var_snp = 2 * maf * (1 - maf)  # n x 1
  var_genetics = t(beta ^ 2) %*% var_snp  # k x 1
  var_e = (1 - h2) * var_genetics / h2  # k x 1
  env_noise = matrix(rnorm(n * k, sd = sqrt(rep(var_e, n))), byrow = TRUE, ncol = k, nrow = n)  # n x k
  pheno = h[[1]] %*% beta + h[[2]] %*% beta + env_noise
  pheno
}

simulate_pheno_single_snp = function(h, beta, h2, maf, null = FALSE) {
  # h: matrix n x p
  # beta: vector p x 1
  # h2: vector (or scalar) p x 1
  # maf: vector p x 1
  n = nrow(h[[1]])
  p = ncol(h[[1]])
  var_snp = 2 * maf * (1 - maf)  # p x 1
  var_genetics = beta ^ 2 * var_snp  # p x 1
  if(null == FALSE) {
    var_e = (1 - h2) * var_genetics / h2  # p x 1
  } else if(null == TRUE) {
    var_e = rep(1, p)  # p x 1
  }
  # print(var_genetics[1])
  # print(var_e[1])
  # print(h2[1])
  env_noise = matrix(rnorm(n * p, sd = sqrt(rep(var_e, n))), byrow = TRUE, ncol = p, nrow = n)  # n x k
  pheno = h[[1]] %*% diag(beta) + h[[2]] %*% diag(beta) + env_noise
  pheno
}

rbeta_from_mean_and_sd = function(n, mean_, sd_) {
  co = mean_ * (1 - mean_) / sd_ ^ 2 - 1
  alpha = mean_ * co
  beta = (1 - mean_) * co
  rbeta(n, shape1 = alpha, shape2 = beta)
}
