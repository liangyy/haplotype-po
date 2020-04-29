# test chr16 data pre-formatted at notebook/03_*

rm(list=ls())

df_father = read.csv('~/Desktop/tmp/haplotype-po/notebook/ukb_eur.father.csv.gzip')
df_mother = read.csv('~/Desktop/tmp/haplotype-po/notebook/ukb_eur.mother.csv.gzip')
df_prs1 = read.csv('~/Desktop/tmp/haplotype-po/notebook/ukb_eur.prs1.csv.gzip')
df_prs2 = read.csv('~/Desktop/tmp/haplotype-po/notebook/ukb_eur.prs2.csv.gzip')

source('code/rlib_em_degenerate.R')
source('code/rlib_em.R')

# rearrange rows 
eid = intersect(df_father$eid, df_mother$eid)
df_father = df_father[match(eid, df_father$eid), ]
df_mother = df_mother[match(eid, df_mother$eid), ]
df_prs1 = df_prs1[match(eid, df_prs1$eid), ]
df_prs2 = df_prs2[match(eid, df_prs1$eid), ]

# rearrange columns 
cols = colnames(df_father)
cols = cols[-length(cols)]
cols = cols[-which(cols == 'Alzheimer.s.disease.dementia')]
df_father = df_father[, cols]
df_mother = df_mother[, cols]
df_prs1 = df_prs1[, cols]
df_prs2 = df_prs2[, cols]

# remove -2
not_minus2 = rowSums(df_father[, -1] == -2) == 0 & rowSums(df_mother[, -1] == -2) == 0
df_father = df_father[not_minus2, ]
df_mother = df_mother[not_minus2, ]
df_prs1 = df_prs1[not_minus2, ]
df_prs2 = df_prs2[not_minus2, ]

# I've implemented the procedure to allow intercept, no need to do centering 
# # center columns
# df_father[, -1] = sweep(df_father[, -1], 2, colMeans(df_father[, -1]), FUN = '-')
# df_mother[, -1] = sweep(df_mother[, -1], 2, colMeans(df_mother[, -1]), FUN = '-')
# 
# # standardize columns
# df_father[, -1] = sweep(df_father[, -1], 2, apply(df_father[, -1], 2, sd), FUN = '/')
# df_mother[, -1] = sweep(df_mother[, -1], 2, apply(df_mother[, -1], 2, sd), FUN = '/')

# run EM
o_deg = em_algorithm_deg(df_father[, -1], df_mother[, -1], df_prs1[, -1], df_prs2[, -1])
o = em_algorithm(df_father[, -1], df_mother[, -1], df_prs1[, -1], df_prs2[, -1])

# Pr(Z) vs number of informative phenotypes
library(dplyr)
library(ggplot2)
n_informative = rowSums(sign(df_father[, -1]) != sign(df_mother[, -1]))
rbind(
  data.frame(n_informative = n_informative, prob_z = o$z_prob_n, type = 'full'),
  data.frame(n_informative = n_informative, prob_z = o_deg$z_prob_n, type = 'degenerate')
) %>% ggplot() + geom_boxplot(aes(x = factor(n_informative), y = abs(prob_z - 0.5), color = type))
table(n_informative)
# hist(o$z_prob_n)

data.frame(n_informative = n_informative, prob_z = o$z_prob_n, prob_z_deg = o_deg$z_prob_n) %>% 
  ggplot() + geom_point(aes(x = prob_z, y = prob_z_deg, color = n_informative)) + coord_equal()

out = data.frame(eid = df_father$eid, prob_z = o$z_prob_n, prob_z_deg = o_deg$z_prob_n)

write.table(out, 'analysis_output/test_em_output.tsv', quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')
