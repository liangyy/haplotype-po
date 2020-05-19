# rm(list = ls())
n = 1000
k = 10
p = 5
h1 = matrix(rnorm(n * k), ncol = k, nrow = n)
h2 = matrix(rnorm(n * k), ncol = k, nrow = n)
beta_f = matrix(rnorm(k * p), ncol = p, nrow = k)
beta_m = matrix(rnorm(k * p), ncol = p, nrow = k)
yf = h1 %*% beta_f + matrix(rnorm(n * p), ncol = p, nrow = n) * 10
ym = h2 %*% beta_m + matrix(rnorm(n * p), ncol = p, nrow = n) * 10
# # 
# library(reticulate)
# np <- import("numpy")
# h1 = np$load('notebook/h1.npy')
# h2 = np$load('notebook/h2.npy')
# yf = np$load('notebook/yf.npy')
# ym = np$load('notebook/ym.npy')

source('code/rlib_em_per_snp.R')
o = em_per_snp(yf, ym, h1[, 1:5], h2[, 1:5])
o$inter
o$beta
o$lld
o$sigma2

hist(o$prob_z)

# o$lld[20:68] + 972676

l1f = NULL
l0f = 0
l1m = NULL
l0m = 0
for(i in 1 : p) {
  yfp1 = .mat_vec_by_row(.mat_vec_by_row(h1, o$beta$f[, i], '*'), o$inter$f[, i], '+')
  yfp2 = h2 * o$beta$f[, i] + o$inter$f[, i]
  ymp1 = .mat_vec_by_row(.mat_vec_by_row(h2, o$beta$m[, i], '*'), o$inter$m[, i], '+')
  ymp2 = h1 * o$beta$m[, i] + o$inter$m[, i]
  rfp1 = .mat_vec_by_col(yfp1, yf[, i], '-'); rfp1_sq = .mat_vec_by_row(rfp1 ^ 2, o$sigma2$f[, i], '/')
  rfp2 = yf[, i] - yfp2; rfp2_sq = rfp2 ^ 2 / o$sigma2$f[, i]
  rmp1 = .mat_vec_by_col(ymp1, ym[, i], '-'); rmp1_sq = .mat_vec_by_row(rmp1 ^ 2, o$sigma2$m[, i], '/')
  rmp2 = ym[, i] - ymp2; rmp2_sq = rmp2 ^ 2 / o$sigma2$m[, i]
  l1f = cbind(l1f, rowSums(rfp1_sq)); l0f = l0f + rfp1_sq
  l1m = cbind(l1m, rowSums(rmp1_sq)); l0m = l0m + rmp2_sq
}
# sum(l1)
# sum(l0)

lf = get_l(yf, h1, h2, o$beta$f, o$inter$f, o$sigma2$f)
lm = get_l(ym, h2, h1, o$beta$m, o$inter$m, o$sigma2$m)
get_lld(lf, lm, 0.5, o$sigma2$f, o$sigma2$m)

l1 = rowSum(l1) + log(0.5)
l0 = rowSum(l0) + log(0.5)
.logsum(-sum(l1), -sum(l0)) - n * sum(log(o$sigma2$f) + log(o$sigma2$m))


plot(o$lld[2:30])

# y = yf[, 2]
# x = h1[, 2]
# 
# mod = lm(y ~ x)
# logLik(mod)
# 
# b = mod$coefficients
# 
# yp = x * b[2] + b[1]
# s2 = summary(mod)$sigma^2
# 
# -  sum((yp - y)^2) / 2 / s2 - 0.5 * length(x) * log(2 * pi * s2)

source('code/rlib_em_otf.R')
o2 = em_algorithm_otf(yf, ym, h1, h2)
# # o$prob_z
# plot(o$lld[5:20])
# hist(o$prob_z)
plot(o2$lld[2:20])
hist(o2$gamma)
hist(o$prob_z)
plot(o$prob_z, o2$gamma)
# 
# source('code/rlib_em_per_snp.R')
# x1 = rnorm(n)
# x2 = rnorm(n)
# y1 = x1 + rnorm(n)
# y2 = x2 + rnorm(n)
# w = runif(n / 2)
# 
# lm_wrap = function(y, x, w) {
#   mod = lm(y ~ x, weights = w)
#   beta = mod$coefficients
#   resid = y - cbind(rep(1, length(x)), x) %*% as.numeric(beta)
#   sigma2 = sum(resid ^ 2 * w) / sum(w)     
#   list(beta = beta, sigma2 = sigma2)
# }
# 
# b = update_beta_and_inter(cbind(x1, x2), cbind(y1, y2), w)
# update_sigma2(cbind(x1, x2), cbind(y1, y2), w, b$beta, b$inter)
# b
# 
# lm_wrap(y1, x1, c(w, 1-w))
# lm_wrap(y2, x2, c(w, 1-w))
# lm_wrap(y2, x1, c(w, 1-w))
