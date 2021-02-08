library(doParallel)

# get list of all cobinding matrices
cobindDir = "/u/scratch/s/seanchea/cobinding/cobind/"
maps <- list.files(cobindDir)
maps <- maps[!(maps %in% c('noBinding.csv'))]
setwd(cobindDir)

# 36 threads is very overkill
cl <- makeCluster(36, type = "PSOCK")
registerDoParallel(cl)
factors <- foreach(matrix = maps, .packages = c("data.table", "pheatmap", "dendsort", "RColorBrewer", "fastcluster")) %dopar% {
    # wrapper for dendsort
    sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

    # remove DE gene names, convert to matrix for pheatmap, transpose
    mat <- fread(matrix, sep = "\t")
    genes <- mat[,1]
    mat <- as.matrix(mat[,-1])
    mat <- t(mat)

    name <- split(matrix, ".")[[1]][1]
    pdf(paste0(cobindDir,name,".pdf"))
    pheatmap(
        mat               = mat,
        cluster_cols      = sort_hclust(hclust.vector(t(mat), method="ward")),
        cluster_rows      = sort_hclust(hclust.vector(mat, method="ward")),
        fontsize          = 12,
        treeheight_row    = 3, 
        treeheight_col    = 3,
        show_colnames     = F,
        fontsize_row      = 4,
        color             = colorRampPalette(brewer.pal(9,"Reds"))(400))
    dev.off()
}
stopCluster(cl)