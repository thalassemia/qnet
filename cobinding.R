library(tidyverse)
library(readr)
library(grid)
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(fastcluster)

intersectDir <- system("echo $SCRATCH/output/intersect/", intern = TRUE)
outDir <- system("echo $SCRATCH/output/cobinding/", intern = TRUE)
genes <- data.frame("rank" = c(1:20, 1:20), 
                    "UD" = c(rep("Down", 20), rep("Up", 20)), 
                    "name" = c("DLX2", "DLX1", "DLX3", "ZNF695", "ATOH8", "MYB", 
                               "NKX3-1", "MYBL2", "GLI1", "CENPA", "SP6", 
                               "ZNF850", "FOXM1", "ZNF492", "E2F1", "GATA2", 
                               "E2F8", "E2F7", "KLF10", "FOSL1", "IRX6", 
                               "PKNOX2", "EGR2", "IRF4", "NR3C2", "THRB", "OSR2",
                               "TSHZ2", "KLF15", "MAF", "GRHL1", "BATF2", "EGR3",
                               "SIX2", "ZNF540", "EBF4", "KLF4", "SATB1", "ZMAT1",
                               "ZNF608"))
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
taskid <- Sys.getenv("SGE_TASK_ID")
set.seed(10)
tfchunks <- split(genes, sample(1:40))
genes <- tfchunks[[taskid]]
rm(tfchunks, taskid)

cobinding <- function(tf) {
  setwd(paste(intersectDir, tf, sep=""))
  intersectBeds <- list.files()
  ca <- NULL

  progress <- 0

  for (bed in intersectBeds) {
    overlap <- read_csv(bed, col_names = FALSE, col_types = "ccccc")
    overlap <- select(overlap, c(4,5))
    #overlap[overlap=="."] <- "0"
    overlap[overlap=="-1.0"] <- "0"
    overlap[,1] <- lapply(overlap[,1], parse_number)
    overlap <- arrange(overlap, overlap[,1], overlap[,2])
    overlap <- distinct(overlap, overlap[,1], .keep_all = TRUE) %>% t()
    row <- shQuote(overlap[1,], type = "cmd")
    row <- paste(row, as.character(unlist(overlap[2,])), sep="=")
    row <- paste(row, collapse=", ")
    temp <- eval(parse(text = paste("data.frame(", row, ")")))
    ca <- rbind(ca, temp)
    intersectBeds[intersectBeds==bed] <- strsplit(bed, "_")[[1]][2]
    progress <- progress + 1
    print(paste(progress/length(intersectBeds)*100, "%", sep=""))
  }
  rownames(ca) <- intersectBeds
  rm(overlap, temp, row, bed, intersectBeds)
  print(paste("Generating heatmap for", tf))
  print("Without heatmap")
  print(gc(full=TRUE))
  pdf(paste(outDir,tf,".png", sep=""))
  pheatmap(
    mat               = as.matrix(ca),
    cluster_cols      = sort_hclust(hclust.vector(t(ca), method="ward")),
    cluster_rows      = sort_hclust(hclust.vector(ca, method="ward")),
    fontsize          = 20,
    treeheight_row    = 0, 
    treeheight_col    = 0,
    show_colnames     = F,
    fontsize_row      = 7,
    main              = paste(tf, " Co-Binding Map (", genes$UD[genes$name==tf], " #", genes$rank[genes$name==tf], ")", sep=""),
    color             = colorRampPalette(brewer.pal(9,"Reds"))(400))
  dev.off()
  print("Heatmap done and with matrix")
  print(gc(full=TRUE))
  rm(ca)
  print("Heatmap done without matrix")
  print(gc(full=TRUE))
}

sigDir <- paste(intersectDir, "signal/", sep="")
qDir <- paste(intersectDir, "qval/", sep="")
sigOut <- paste(outDir, "signal/", sep="")
qOut <- paste(outDir, "qval/", sep="")
intersectDir <- sigDir
outDir <- sigOut
tfs <- list.files(intersectDir)
tf <- tfs[tfs %in% genes$name]
rm(tfs)
stopifnot(!is.null(tf))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
cobinding(tf)
intersectDir <- qDir
outDir <- qOut
tfs <- list.files(intersectDir)
tf <- tfs[tfs %in% genes$name]
rm(tfs)
stopifnot(!is.null(tf))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
cobinding(tf)