library(tidyverse)

deFile <- ""
outDir <- ""
targetsDir <- ""

de <- read_delim(deFile, "\t") %>% rename(GeneSymbol = gene_name) %>% filter(padj <= 0.05)

setwd(targetsDir)
TFs <- list.files()

for(TF in TFs) {
  setwd(targetsDir)
  targets <- read_delim(TF, "\t")
  targets <- inner_join(targets, de, by="GeneSymbol")
  setwd(outDir)
  write.table(targets, TF, sep=",")
}