library(tidyverse)
library(fastcluster)
library(pheatmap)
library(RColorBrewer)

detargets = "/u/scratch/s/seanchea/detargets"
outDir = "/u/home/s/seanchea/"
setwd(detargets)

tfs <- list.files(detargets)
matrix <- read_csv(tfs[1], col_types=cols()) %>% select(c(5,7,10))

for (tf in tfs[-1]) {
  temp <- read_csv(tf, col_types=cols()) %>% select(c(5,7,10))
  matrix <- full_join(matrix, temp, by="GeneSymbol")
}

score <- NULL
targets <- matrix[2]
matrix <- matrix[-2]
matrix[is.na(matrix)] <- 0


for (i in 1:length(tfs)) {
  score <- rbind(score, t(matrix[,1]*matrix[,2]))
}

rownames(score) <- unlist(strsplit(tfs,".csv"))
colnames(score) <- targets

png(paste(outDir,"enrich.png", sep=""), width=1200, height=1200, res=300)
pheatmap(
  mat               = matrix,
  cluster_cols      = sort_hclust(hclust.vector(t(matrix), method="ward")),
  cluster_rows      = sort_hclust(hclust.vector(matrix, method="ward")),
  fontsize          = 20,
  treeheight_row    = 0, 
  treeheight_col    = 0,
  show_colnames     = F,
  fontsize_row      = 7,
  main              = paste(tf, " Co-Binding Map (", genes$UD[genes$name==tf], " #", genes$rank[genes$name==tf], ")", sep=""),
  color             = colorRampPalette(brewer.pal(9,"Reds"))(400))
dev.off()