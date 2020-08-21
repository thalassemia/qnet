library(tidyverse)
library(readr)
library(grid)
library(pheatmap)
library(dendsort)
library(RColorBrewer)

dir <- "/home/sean/TFChIP/intersect/"
outDir <- "/home/sean/TFChIP/cobinding/"
genes <- data.frame("rank" = c(1:20, 1:20), 
                    "UD" = c(rep("Down", 20), rep("Up", 20)), 
                    "name" = c("DLX2", "DLX1", "DLX3", "ZNF695", "ATOH8", "MYB", 
                               "NKX3-1", "MYBL2", "GL1", "CENPA", "SP6", 
                               "ZNF850", "FOXM1", "ZNF492", "E2F1", "GATA2", 
                               "E2F8", "E2F7", "KLF10", "FOSL1", "IRX6", 
                               "PKNOX2", "EGR2", "IRF4", "NR3C2", "THRB", "OSR2",
                               "TSHZ2", "KLF15", "MAF", "GRHL1", "BATF2", "EGR3",
                               "SIX2", "ZNF540", "EBF4", "KLF4", "SATB1", "ZMAT1",
                               "ZNF608"))

setwd(dir)
tfs <- list.files()
tfs <- tfs[tfs %in% genes$name]
# tfs <- tfs[!(tfs %in% c("FOSL1", "EGR3"))]

progress <- 0

for (tf in tfs) {
  tf <- "EGR3"
  setwd(file.path(dir, tf))
  intersectBeds <- list.files()
  ca <- NULL
  
  indivProg <- 0
  
  for (bed in intersectBeds) {
    overlap <- read_csv(bed, col_names = FALSE, , col_types = cols())
    overlap <- select(overlap, c(4,ncol(overlap)))
    overlap[overlap=="."] <- "0"
    overlap[,1] <- lapply(overlap[,1], as.character)
    overlap[,1] <- lapply(overlap[,1], parse_number)
    overlap <- distinct(overlap, X4, .keep_all = TRUE)
    overlap <- arrange(overlap, select(overlap,1), select(overlap,2)) %>% t()
    overlap[1,] <- as.numeric(overlap[1,])
    row <- as.vector(outer("peak", overlap[1,], paste, sep=""))
    row <- shQuote(row, type = "cmd")
    row <- paste(row, as.character(unlist(overlap[2,])), sep="=")
    row <- paste(row, collapse=", ")
    temp <- tryCatch(
      {
        eval(parse(text = paste("data.frame(", row, ")")))
      }, 
      error = function(e) {
       return(-1)
      }
    )
    if (temp == -1) {
      intersectBeds = intersectBeds[intersectBeds!=bed]
      indivProg <- indivProg - 1
    }
    else {
      ca <- rbind(ca, temp)
    }
    intersectBeds[intersectBeds==bed] <- strsplit(bed, "_")[[1]][2]
    indivProg <- indivProg + 1
    print(paste("Total", progress/length(tfs)*100, "%; ", tf, indivProg/length(intersectBeds)*100, "%"))
  }
  rownames(ca) <- intersectBeds
  progress <- progress + 1
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  mat_cluster_cols <- sort_hclust(hclust(dist(t(ca))))
  mat_cluster_rows <- sort_hclust(hclust(dist(ca)))
  png(paste(outDir,tf,".jpg", sep=""), width=1000, height=800)
  pheatmap(
    mat               = as.matrix(ca),
    cluster_cols      = mat_cluster_cols,
    cluster_rows      = mat_cluster_rows,
    fontsize          = 20,
    treeheight_row    = 0, 
    treeheight_col    = 0,
    show_colnames     = F,
    fontsize_row      = 7,
    main              = paste(tf, " Co-Binding Map (", genes$UD[genes$name==tf], " #", genes$rank[genes$name==tf], ")", sep=""),
    color             = colorRampPalette(brewer.pal(9,"Reds"))(400))
  dev.off()
  break
}