library(tidyverse)
library(fastcluster)
library(pheatmap)
library(RColorBrewer)

detargets = "/u/scratch/s/seanchea/detargets"
setwd(detargets)

tfs <- list.files(detargets)
matrix <- read_csv(tfs[1], col_types=cols()) %>% select(c(5,7,10,14))

for (tf in tfs[-1]) {
  temp <- read_csv(tf, col_types=cols()) %>% select(c(5,7,10,14))
  matrix <- full_join(matrix, temp, by="GeneSymbol")
}

tfs <- unlist(strsplit(tfs,".csv"))