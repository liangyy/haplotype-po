library(dplyr)
library(data.table)
options(datatable.fread.datatable = FALSE)
library(ggplot2)

# look at GWAS summary statistics from GWAX paper 
# my best guess is that they are using PLINK format so that A1 is effect allele
# reference: https://www.biostars.org/p/90264/
gwas_gwax = fread(cmd = 'zcat < ~/Downloads/AD.meta.assoc.gz', header = TRUE, sep = ' ') 
sum(gwas_gwax$P_META < 1e-50)
gwas_gwax_chr16 = gwas_gwax %>% filter(CHR == 16) 

# IGAP http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php
gwas_igap = fread('~/Downloads/IGAP_summary_statistics/IGAP_stage_1.txt', header = TRUE, sep = '\t') 
gwas_igap_chr16 = gwas_igap %>% filter(Chromosome == 16) %>% mutate(Pvalue = as.numeric(Pvalue))


# gwas_gwax_chr16 %>% ggplot() + geom_point(aes(x = BP, y = -log(P_META)))
snp_map = fread(cmd = 'zcat < ~/Desktop/tmp/haplotype-po/snp_map_for_neale_lab_gwas.with_sign.tsv.gz', header = TRUE, sep = '\t')


# join the two
allele_checker = function(allele_id, a1, a2) {
  # assume a1 is effect allele
  str_p = paste0(a2, ',', a1)
  str_n = paste0(a1, ',', a2)
  o = rep(NA, length(allele_id))
  o[allele_id == str_p] = '+'
  o[allele_id == str_n] = '-'
  return(o)
}
merge = gwas_igap_chr16 %>% inner_join(snp_map, by = c('Position' = 'pos'))
merge = merge %>% mutate(direction_of_match = allele_checker(allele_ids, Effect_allele, Non_Effect_allele))
table(merge$direction_of_match)

hist(merge$Pvalue)
source('https://gist.githubusercontent.com/liangyy/605437935751bb77b1739666b18517bf/raw/85b1e83d1cbf7c8b5671030e927000342ed16f97/my_qqplot_by_group.R')
my_qqplot_by_group(gwas_igap_chr16$Pvalue, '.')
