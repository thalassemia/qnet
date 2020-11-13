database <- read.csv("Database.csv")
names(database)[names(database) == "HGNC.symbol"] <- "gene_name"
SSvsP <- read.csv("SSvsP.csv")
combined <- inner_join(database, SSvsP, by = "gene_name")
combined <- filter(combined, padj <= 0.05 & str_detect(combined$Is.TF., "Yes") & (log2FoldChange >=1 | log2FoldChange <=-1)) %>% arrange(log2FoldChange)
write.csv(combined, "SSvsP_SigTFs.csv")