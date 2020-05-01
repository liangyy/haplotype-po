# x = sort(runif(100))
# y = c(rep(0, 50), rep(1, 50))
# 
# mod = glm(cbind(y, 1 - y) ~ 1 + x, family = 'binomial')

library(dplyr)
library(ggplot2)
library(data.table)

library(reticulate)
np <- import("numpy")
check_sum = np$load('/Users/yanyul/Desktop/tmp/haplotype-po/test_sanity_sum.npy')
check_avg = np$load('/Users/yanyul/Desktop/tmp/haplotype-po/test_sanity_avg.npy')
check_imp = np$load('/Users/yanyul/Desktop/tmp/haplotype-po/test_sanity_imp.npy')
dim(check_sum)
dim(check_avg)
dim(check_imp)

gwas_father = fread('~/Downloads/2_UKB_AD_paternal_summary_output_June2019.txt', sep = ' ', data.table = FALSE)
snp_map = fread('zcat < ~/Desktop/tmp/haplotype-po/snp_map_for_neale_lab_gwas.with_sign.tsv.gz', header = TRUE, sep = '\t')
snp_map_16 = snp_map %>% filter(chrom == 16)

gwas_father_extracted = left_join(snp_map_16, gwas_father %>% filter(CHR == 16), by = c('pos' = 'BP'))

to_df = function(mat) {
  df = as.data.frame(mat)
  colnames(df) = c('beta', 'se', 'convergence')
  return(df)
}
df_trait1 = inner_join(
  to_df(t(check_sum[, 1, ])) %>% mutate(snpid = 1 : nrow(.)),
  to_df(check_avg[1, 1, , ]) %>% mutate(snpid = 1 : nrow(.)),
  by = 'snpid',
  suffix = c('_sum', '_avg')
)
df_trait2 = inner_join(
  to_df(t(check_sum[, 2, ])) %>% mutate(snpid = 1 : nrow(.)),
  to_df(check_avg[1, 2, , ]) %>% mutate(snpid = 1 : nrow(.)),
  by = 'snpid',
  suffix = c('_sum', '_avg')
)
df_trait2_cont = inner_join(
  to_df(check_imp[1, 2, , ]) %>% mutate(snpid = 1 : nrow(.)),
  to_df(check_imp[2, 2, , ]) %>% mutate(snpid = 1 : nrow(.)),
  by = 'snpid',
  suffix = c('_imp', '_flip')
)
df_trait2 = inner_join(df_trait2, df_trait2_cont, by = 'snpid')
gwas_father_extracted = cbind(
  gwas_father_extracted, df_trait2
)

# convergence
gwas_father_extracted %>% group_by(convergence_sum, convergence_avg, convergence_imp, convergence_flip) %>% summarize(n())
gwas_father_extracted_converged = gwas_father_extracted %>% 
  filter(convergence_sum == 1, convergence_avg == 1, convergence_imp == 1, convergence_flip == 1) %>% 
  filter(!is.na(se_sum), !is.na(se_avg), se_sum > 0, se_avg > 0) %>% 
  filter(!is.na(se_imp), !is.na(se_flip), se_imp > 0, se_flip > 0)

# take common
gwas_father_extracted_converged = gwas_father_extracted_converged %>% filter(!is.na(CHR))

# fix sign
signed_beta = function(beta, allele_ids, a1, a2) {
  ref = unlist(lapply(strsplit(allele_ids, ','), function(x) {x[1]}))
  alt = unlist(lapply(strsplit(allele_ids, ','), function(x) {x[2]}))
  obeta = rep(NA, length(beta))
  obeta[ref == a2 & alt == a1] = beta[ref == a2 & alt == a1]
  obeta[ref == a1 & alt == a2] = -beta[ref == a1 & alt == a2]
  return(obeta)
}
gwas_father_extracted_converged = gwas_father_extracted_converged %>% mutate(corrected_beta = signed_beta(Beta, allele_ids, A1, A2))

gwas_father_extracted_converged %>% ggplot() + 
  geom_point(aes(x = beta_avg, y = beta_sum, color = 'beta_sum'), alpha = 0.2) +
  geom_point(aes(x = beta_avg, y = corrected_beta, color = 'beta_paper'), alpha = 0.2) +
  geom_point(aes(x = beta_avg, y = beta_imp, color = 'beta_imp'), alpha = 0.2) + 
  geom_point(aes(x = beta_avg, y = beta_flip, color = 'beta_flip'), alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0) + coord_equal() + ylab('beta (sum or paper)')

gwas_reshape = inner_join(
  gwas_father_extracted_converged %>% 
    select(pos, corrected_beta, beta_avg, beta_sum, beta_imp, beta_flip) %>% 
    rename(paper = corrected_beta, avg = beta_avg, sum = beta_sum, imp = beta_imp, flip = beta_flip) %>%
    reshape2::melt(id.vars = c('pos')),
  gwas_father_extracted_converged %>% 
    select(pos, SE, se_avg, se_sum, se_imp, se_flip) %>% 
    rename(paper = SE, avg = se_avg, sum = se_sum, imp = se_imp, flip = se_flip) %>% 
    reshape2::melt(id.vars = c('pos')),
  by = c('pos', 'variable'),
  suffix = c('_beta', '_se')
) %>% mutate(pval = exp(pnorm(-abs(value_beta / value_se), log.p = T)) * 2) 

gwas_reshape = gwas_reshape %>% group_by(variable) %>% mutate(pexp = rank(pval) / (n() + 1))
gwas_reshape %>% ggplot() + 
  geom_point(aes(x = -log(pexp), y = -log(pval))) + 
  facet_wrap(~variable, ncol = 3) + geom_abline(slope = 1, intercept = 0) + coord_equal()

gwas_reshape_reshape = gwas_reshape %>% reshape2::dcast(pos ~ variable, value.var = 'pval') 
gwas_reshape_reshape %>% ggplot() +
  geom_point(aes(x = -log(paper), y = -log(avg)), alpha = 0.2) + geom_abline(slope = 1, intercept = 0, color = 'red')

gwas_reshape_reshape %>% ggplot() +
  geom_point(aes(x = -log(avg), y = -log(imp)), alpha = 0.2) + geom_abline(slope = 1, intercept = 0, color = 'red')

gwas_reshape_reshape %>% ggplot() +
  geom_point(aes(x = -log(avg), y = -log(flip)), alpha = 0.2) + geom_abline(slope = 1, intercept = 0, color = 'red')



tops = gwas_reshape_reshape %>% arrange(desc(-avg)) %>% head(2)
gwas_father_extracted_converged %>% filter(pos %in% tops$pos)
