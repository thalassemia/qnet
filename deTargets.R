library(tidyverse)

deFile <- "/home/sean/data/de/deseq_SSvsP_gencodev29_allgenes_021120.txt"
outDir <- "/home/sean/data/detargets/"
targetsDir <- "/home/sean/data/targets/"

de <- read_delim(deFile, "\t") %>% rename(GeneSymbol = gene_name) %>% filter(padj <= 0.05)

setwd(targetsDir)
TFs <- list.files()

for(TF in TFs) {
  setwd(targetsDir)
  targets <- read_delim(TF, "\t", col_types = cols())
  targets <- inner_join(targets, de, by="GeneSymbol")
  targets <- arrange(targets, desc(log2FoldChange))
  setwd(outDir)
  write.table(targets, TF, sep=",", row.names = F)
}