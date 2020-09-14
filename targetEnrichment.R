library(tidyverse)
library(fastcluster)
library(pheatmap)
library(RColorBrewer)
library(dendsort)
library(parallel)

detargets = system("echo $SCRATCH/output/detargets", intern = TRUE)
outDir = system("echo $SCRATCH/output/enrich/", intern = TRUE)
de <- system("echo $SCRATCH/data/deTF.csv", intern = TRUE)
t <- seq(0,3,length.out=7)
dir.create(outDir, showWarnings=FALSE)
setwd(detargets)
de <- read_csv(de, col_types=cols()) %>% select(c(gene_name, log2FoldChange))
de <- transmute(de, tfFiles = paste0(de$gene_name, rep(".csv", nrow(de))), up = log2FoldChange > 0)
allfactors <- list.files(detargets)
de <- filter(de, tfFiles %in% allfactors)
direction <- function(thresh) {
  uptfs <- filter(de, up) %>% select(tfFiles)
  downtfs <- filter(de, !up) %>% select(tfFiles)
  enrich(unlist(uptfs), "Up", thresh)
  print(paste("Done with Up", thresh))
  enrich(unlist(downtfs), "Down", thresh)
  print(paste("Done with Down", thresh))
}
enrich <- function(tfs, direc, threshold) {
  matrix <- read_csv(tfs[1], col_types=cols()) %>% select(c(Score, GeneSymbol, log2FoldChange))
  matrix <- filter(matrix, Score>=threshold)
  matrix$factor <- rep(strsplit(tfs[1],".csv"), nrow(matrix))
  
  for (tf in tfs[-1]) {
    temp <- read_csv(tf, col_types=cols()) %>% select(c(Score, GeneSymbol, log2FoldChange))
    temp$factor <- rep(strsplit(tf,".csv"), nrow(temp))
    temp <- filter(temp, Score>=threshold)
    matrix <- rbind(matrix, temp)
  }
  
  targets <- matrix[!duplicated(matrix[,2]),2]
  
  scores <- matrix(0, length(tfs), nrow(targets))
  rownames(scores) <- unlist(strsplit(tfs,".csv"))
  colnames(scores) <- targets$GeneSymbol
  for (row in 1:nrow(matrix)) {
    scores[unlist(matrix[row,4]), unlist(matrix[row,2])] = unlist(matrix[row,3])
  }
  
  count <- data.frame(up=rowSums((scores>0)*1), down=rowSums((scores<0)*1))
  rownames(count) <- unlist(strsplit(tfs,".csv"))
  print(count)
  print(count$up/count$down)
  count[["Ratio of Up to Down"]] <- count$up/count$down
  temp <- scores
  temp[temp==0] <- NA
  count$Mean <- rowMeans(temp,na.rm=TRUE)
  rm(temp)
  probabilities <- c(0.05, 0.95)
  percentiles <- apply(scores,1,function(x) quantile(x[x!=0],probs=probabilities)) %>% t()
  count <- cbind(count, percentiles)
  count <- arrange(count, desc(up))
  count$"Up Target Rank" <- 1:nrow(count)
  count <- arrange(count,desc(down))
  count$"Down Target Rank" <- 1:nrow(count)
  count <- arrange(count, desc("Ratio of Up to Down"))
  count$"Ratio of Up to Down Rank" <- 1:nrow(count)
  count <- arrange(count, desc(Mean))
  count$"Mean Rank" <- 1:nrow(count)
  for (i in colnames(percentiles)) {
    count <- arrange(count, desc(count[[i]]))
    count[[paste(i,"Rank")]] <- 1:nrow(count)
  }
  write.table(count, file=paste(outDir, direc, "count", threshold,".csv",sep=""), sep=",", col.names=NA)
  write.table(scores, file=paste(outDir, direc, "enrich", threshold,".csv", sep=""), sep=",", col.names=NA)
  
  paletteLength <- 1000
  myBreaks <- c(seq(min(scores), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(scores)/paletteLength, max(scores), length.out=floor(paletteLength/2)))
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  clustcols <- sort_hclust(hclust.vector(t(scores), method="ward"))
  clustrows <- sort_hclust(hclust.vector(scores, method="ward"))
  png(paste(outDir,direc,"enrich",threshold,".png", sep=""), width=5, height=5, units="in", res=1000, pointsize=12)
  pheatmap(
    mat               = scores,
    cluster_cols      = clustcols,
    cluster_rows      = clustrows,
    fontsize          = 20,
    treeheight_row    = 0, 
    treeheight_col    = 0,
    show_colnames     = F,
    fontsize_row      = 7,
    main              = paste("Enrichment for Putative Targets of ", direc, "regulated TFs", sep=""),
    color             = colorRampPalette(brewer.pal(11,"Spectral"))(paletteLength),
    breaks            = myBreaks)
  dev.off()
}
mclapply(t, direction)