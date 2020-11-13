library(data.table)
library(tidyverse)
library(doParallel)
library(grid)
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(fastcluster)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
cobindDir = "/u/scratch/s/seanchea/tfchip/Lcobind/"

ca <- fread(file.path(cobindDir, "lisa_train.csv"), sep = "\t")
genes <- ca[,1]
#ca <- as.data.frame(t(ca[,-1]))
#colnames(ca) <- genes[[1]]
ca <- as.data.frame(ca[-1,-1])
rownames(ca) <- genes[[1]][-1]

pdf(paste0(cobindDir,"lisa_train.pdf"))
pheatmap(
    mat               = as.matrix(ca),
    cluster_cols      = sort_hclust(hclust.vector(t(ca), method="ward")),
    cluster_rows      = sort_hclust(hclust.vector(ca, method="ward")),
    fontsize          = 12,
    treeheight_row    = 0, 
    treeheight_col    = 0,
    show_colnames     = F,
    fontsize_row      = 4,
    main              = "Lisa TF Co-Binding Map",
    color             = colorRampPalette(brewer.pal(9,"Reds"))(400))
dev.off()