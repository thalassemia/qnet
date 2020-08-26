library(tidyverse)
library(fastcluster)
library(pheatmap)
library(RColorBrewer)
library(dendsort)

detargets = system("echo $SCRATCH/output/detargets", intern = TRUE)
outDir = system("echo $SCRATCH/output/", intern = TRUE)
dir.create(outDir, warnings=FALSE)
setwd(detargets)

tfs <- list.files(detargets)
matrix <- read_csv(tfs[1], col_types=cols()) %>% select(c(Score, GeneSymbol, log2FoldChange))
matrix$factor <- rep(strsplit(tfs[1],".csv"), nrow(matrix))

for (tf in tfs[-1]) {
  temp <- read_csv(tf, col_types=cols()) %>% select(c(Score, GeneSymbol, log2FoldChange))
  temp$factor <- rep(strsplit(tf,".csv"), nrow(temp))
  matrix <- rbind(matrix, temp)
}
print(gc())

matrix$enrich <- matrix[,3]
targets <- matrix[!duplicated(matrix[,2]),2]

scores <- matrix(0, length(tfs), nrow(targets))
rownames(scores) <- unlist(strsplit(tfs,".csv"))
colnames(scores) <- targets$GeneSymbol

for (row in 1:nrow(matrix)) {
  scores[unlist(matrix[row,4]), unlist(matrix[row,2])] = unlist(matrix[row,5])
}


sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
clustcols <- sort_hclust(hclust.vector(t(scores), method="ward"))
clustrows <- sort_hclust(hclust.vector(scores, method="ward"))
png(paste(outDir,"enrich.png", sep=""), width=1200, height=1200)
pheatmap(
  mat               = scores,
  cluster_cols      = clustcols,
  cluster_rows      = clustrows,
  fontsize          = 20,
  treeheight_row    = 0, 
  treeheight_col    = 0,
  show_colnames     = F,
  fontsize_row      = 7,
  main              = "Enrichment for Putative TF Targets",
  color             = colorRampPalette(brewer.pal(11,"Spectral"))(1000))
dev.off()