library(data.table)
library(tidyverse)
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(fastcluster)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
cobindDir = "/u/scratch/s/seanchea/hmchip/cobind/"

ca <- fread(file.path(cobindDir, "atac.csv"), sep = "\t")
ca <- t(ca[,-1])

pdf(paste0(cobindDir,"atac.pdf"))
pheatmap(
    mat               = as.matrix(ca),
    cluster_cols      = hclust.vector(t(ca), method="ward"),
    cluster_rows      = hclust.vector(ca, method="ward"),
    #cluster_cols      = sort_hclust(hclust.vector(t(ca), method="ward")),
    #cluster_rows      = sort_hclust(hclust.vector(ca, method="ward")),
    fontsize          = 12,
    treeheight_row    = 0, 
    treeheight_col    = 0,
    show_colnames     = F,
    fontsize_row      = 4,
    main              = "DA Region Co-Binding Map",
    color             = colorRampPalette(brewer.pal(9,"Reds"))(400))
dev.off()