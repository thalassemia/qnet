library(tidyverse)

deFile <- system("echo $SCRATCH/data/deseq_SSvsP_gencodev29_allgenes_021120.txt", intern = TRUE)
targetsDir <- system("echo $SCRATCH/output/targets/", intern = TRUE)
outDir <- system("echo $SCRATCH/output/detargets/", intern = TRUE)
dir.create(outDir, showWarnings = FALSE)

de <- read_delim(deFile, "\t", col_types=cols()) %>% rename(GeneSymbol = gene_name) %>% filter(padj <= 0.05)

setwd(targetsDir)
TFs <- list.files()

for(TF in TFs) {
  setwd(targetsDir)
  targets <- read_delim(TF, "\t", col_types = cols())
  targets <- inner_join(targets, de, by="GeneSymbol")
  setwd(outDir)
  write.table(targets, TF, sep=",", row.names = F)
}