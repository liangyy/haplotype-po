# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("huex10sttranscriptcluster.db")

# adapted from https://github.com/hwheeler01/trans-px/blob/master/02_annotate_FHSexp.R
args <- commandArgs(trailingOnly = TRUE)
print(args)

annCols = c("PROBEID", "ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")

library(huex10sttranscriptcluster.db)
library(data.table)
# library(dplyr)
idList = keys(huex10sttranscriptcluster.db)

idTable = select(huex10sttranscriptcluster.db,
                 keys=idList, keytype="PROBEID", columns=annCols)

geneIDs = fread(paste0('cat ', args[1], ' | grep -v "^#%"| cut -f 1'), header = TRUE, data.table = FALSE)
geneIDs$probeset_id = as.character(as.numeric(geneIDs$probeset_id))
# geneMerge = merge(geneIDs, idTable, by.x = "probeset_id", by.y = "PROBEID", all.x = TRUE)
geneMerge = dplyr::left_join(geneIDs, idTable, by = c('probeset_id' = 'PROBEID'))

write.table(geneMerge, args[2], 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
