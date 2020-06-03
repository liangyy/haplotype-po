# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("huex10sttranscriptcluster.db")

# adapted from https://github.com/hwheeler01/trans-px/blob/master/02_annotate_FHSexp.R
args <- commandArgs()


annCols = c("PROBEID", "ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")

library(huex10sttranscriptcluster.db)
library(data.table)

idList = keys(huex10sttranscriptcluster.db)

idTable = select(huex10sttranscriptcluster.db,
                 keys=idList, keytype="PROBEID", columns=annCols)

geneIDs = fread(paste0('cat ', args[1], ' | grep -v "^#%"| cut -f 1'), header = TRUE)

geneMerge = merge(geneIDs, idTable, by.x="probeset_id", by.y="PROBEID", all.x=TRUE)

write.table(geneMerge, args[2], 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
