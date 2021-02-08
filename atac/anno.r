library(tidyverse)
#library(ChIPseeker)
#library(org.Hs.eg.db)
#library(GenomicFeatures)

outDir <- system("echo $SCRATCH/atac/corr/", intern = TRUE)
data <- system("echo $SCRATCH/data/SS_P_diff_bound_sites.csv", intern = TRUE)
gencodeDir <- system("echo $SCRATCH/data/", intern = TRUE)
setwd(gencodeDir)
#data <- readPeakFile(data, sep = ",")
#txdb <- makeTxDbFromGFF("gencode.v29.annotation.gtf.gz", format = "gtf")
#suppressWarnings(dataAnno <- annotatePeak(data, TxDb = txdb, annoDb = "org.Hs.eg.db"))
#dataAnno <- as.data.frame(dataAnno)
setwd(outDir)
#write.csv(dataAnno, "atacAnno.csv", row.names = F)

dataAnno <- read_csv("atacAnno.csv")
dataAnno <- filter(dataAnno, abs(distanceToTSS) <= 3000)
#dePath <- system("echo $SCRATCH/data/bsf.csv", intern = T)
#de <- read_delim(dePath, delim = "\t", col_names = F)
dePath <- system("echo $SCRATCH/data/deseq_SSvsP_gencodev29_allgenes_021120.txt", intern = T)
de <- read_delim(dePath, delim = "\t")
de <- filter(de, padj < 0.05)
de <- dplyr::rename(de, SYMBOL = gene_name)
de <- dplyr::rename(de, rna_de_l2fc = log2FoldChange)
de <- inner_join(de, dataAnno, "SYMBOL")
write.csv(de, "DE.csv", row.names = F)

ggplot(de, aes(x=rna_de_l2fc, y=Fold)) + geom_point() + geom_smooth(method=lm) + labs(x="RNA-Seq Fold Change", y="ATAC-Seq Fold Change")
ggsave("DE_plot.png")

print(shapiro.test(de$rna_de_l2fc))
print(shapiro.test(de$Fold))

pearson <- cor.test(de$rna_de_l2fc, de$Fold, method = "pearson")
correlation <- data.frame(stat = c(pearson$statistic),
df = c(pearson$parameter),
p = c(pearson$p.value),
score = c(pearson$estimate), n = c(nrow(de)))
print(nrow(de))
print(correlation)
rownames(correlation) <- c("Pearson")
write.csv(correlation, "DE_corr.csv")