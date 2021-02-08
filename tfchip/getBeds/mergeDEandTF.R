library(tidyverse)

# both TF list and deseq output are in this directory
setwd(system("echo $SCRATCH/data/", intern = TRUE))

database <- read_csv("Database.csv")
SSvsP <- read_delim("deseq_SSvsP_gencodev29_allgenes_021120.txt", "\t")
# remove version info from Ensembl IDs to maintain intercompatibility
SSvsP$gene_id <- sapply(SSvsP, strsplit(split='.')[[1]][1])
SSvsP <- rename(SSvsP, "Ensembl ID" = "gene_id")
# find DE TFs by merging on Ensembl gene IDs
combined <- inner_join(database, SSvsP, by = "Ensembl ID")
# return file with significant (by p-adj AND log2FoldChange) DE TFs
combined <- filter(combined, padj<=0.05 & str_detect(combined$"Is TF?", "Yes") & abs(log2FoldChange)>=1) %>% arrange(log2FoldChange)
write.csv(combined, "deTF.csv")