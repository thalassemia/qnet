library(data.table)
library(tidyverse)
library(doParallel)
library(grid)
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(fastcluster)

intDir <- system("echo $SCRATCH/tfchip/intersectAnnoSALL2/qval/", intern = TRUE)
oDir <- system("echo $SCRATCH/tfchip/cobindNew/qval/", intern = TRUE)
deData <- system("echo $SCRATCH/data/deseq_SSvsP_gencodev29_allgenes_021120.txt", intern = TRUE)
genes <- data.frame("rank" = c(1:20, 1:20, 0), 
                "UD" = c(rep("Down", 20), rep("Up", 21)), 
                "name" = c("DLX2", "DLX1", "DLX3", "ZNF695", "ATOH8", "MYB", 
                            "NKX3-1", "MYBL2", "GLI1", "CENPA", "SP6", 
                            "ZNF850", "FOXM1", "ZNF492", "E2F1", "GATA2", 
                            "E2F8", "E2F7", "KLF10", "FOSL1", "IRX6", 
                            "PKNOX2", "EGR2", "IRF4", "NR3C2", "THRB", "OSR2",
                            "TSHZ2", "KLF15", "MAF", "GRHL1", "BATF2", "EGR3",
                            "SIX2", "ZNF540", "EBF4", "KLF4", "SATB1", "ZMAT1",
                            "ZNF608", "SALL2"))

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
de <- read_delim(deData, "\t", col_types = cols())

cobinding <- function(tf, intersectDir, outDir) {    
    setwd(paste0(intersectDir, tf))
    files <- list.files()
    names <- c()
    ffranks <- c()
    pfranks <- c()
    for (i in files) {
        split_name <- strsplit(i, "_")[[1]]
        ffranks <- append(ffranks, split_name[1])
        names <- append(names, split_name[2])
        pfranks <- append(pfranks, strsplit(split_name[3], ".csv")[[1]][1]) 
    }
    intersectBeds <- data.frame("file" = files, "name" = names, "ffrank" = ffranks, "pfrank" = pfranks)
    intersectBeds <- intersectBeds %>% arrange(ffrank, pfrank)
    split <- filter(genes, name == tf)
    #intersectBeds <- intersectBeds[temp %in% filter(genes, UD==split$UD[1])$name]
    # Comment out line above and uncomment lines below to include all signicantly enriched/depleted TF's, not just those in the top 20
    if (split$UD == "Up") {
        intersectBeds <- intersectBeds %>% filter(name %in% filter(de, log2FoldChange > 0, padj <= 0.05)$gene_name)
    }
    else {
        intersectBeds <- intersectBeds %>% filter(name %in% filter(de, log2FoldChange < 0, padj <= 0.05)$gene_name)
    }
    ca <- NULL

    progress <- 0
    bestFile <- c()
    partner_factors <- c()
    for (cf in unique(intersectBeds$name)) {
        cf_files <- filter(intersectBeds, name == cf)
        bestFile <- append(bestFile, cf_files$file[1])
        partner_factors <- append(partner_factors, cf)
    }
    for (bed in bestFile) {
        overlap <- fread(bed, sep=",")
        overlap <- select(overlap, c(4,9))
        overlap[overlap=="."] <- "0"
        overlap[,1] <- lapply(overlap[,1], parse_number) %>% lapply(as.numeric)
        overlap <- arrange(overlap, overlap[,1], overlap[,2])
        overlap <- distinct(overlap, overlap[,1], .keep_all = TRUE) %>% t()
        row <- paste0("peak", trimws(overlap[1,]))
        colnames(overlap) <- row
        overlap <- overlap[2,]
        ca <- rbind(ca, overlap)
        progress <- progress + 1
        print(paste("Loading beds for", tf, basename(outDir), progress/length(bestFile)*100, "%"))
    }
    if (!is.null(ca)) {
        rownames(ca) <- partner_factors
        print(paste0(outDir, tf, '_train.csv'))
        write.csv(ca, paste0(outDir, tf, '_train.csv'))
        #pdf(paste(outDir,tf,".pdf", sep=""))
        #pheatmap(
        #    mat               = as.matrix(ca),
        #    cluster_cols      = sort_hclust(hclust.vector(t(ca), method="ward")),
        #    cluster_rows      = sort_hclust(hclust.vector(ca, method="ward")),
        #    fontsize          = 12,
        #    treeheight_row    = 0, 
        #    treeheight_col    = 0,
        #    show_colnames     = F,
        #    fontsize_row      = 10,
        #    main              = paste(tf, " Co-Binding Map (", genes$UD[genes$name==tf], " #", genes$rank[genes$name==tf], ")", sep=""),
        #    color             = colorRampPalette(brewer.pal(9,"Reds"))(400))
        #dev.off()
        rm(ca)
    }

    ca <- NULL

    progress <- 0
    altFile <- c()
    partner_factors <- c()

    for (cf in unique(intersectBeds$name)) {
        cf_files <- filter(intersectBeds, name == cf)
        if (length(cf_files$file) > 1) {
            altFile <- append(altFile, cf_files$file[2])
        }
        else {
            altFile <- append(altFile, cf_files$file[1])
        }
        partner_factors <- append(partner_factors, cf)
    }
    for (bed in altFile) {
        overlap <- fread(bed, sep=",")
        overlap <- select(overlap, c(4,9))
        overlap[overlap=="."] <- "0"
        overlap[,1] <- lapply(overlap[,1], parse_number) %>% lapply(as.numeric)
        overlap <- arrange(overlap, overlap[,1], overlap[,2])
        overlap <- distinct(overlap, overlap[,1], .keep_all = TRUE) %>% t()
        row <- paste0("peak", trimws(overlap[1,]))
        colnames(overlap) <- row
        overlap <- overlap[2,]
        ca <- rbind(ca, overlap)
        progress <- progress + 1
        print(paste("Loading beds for", tf, basename(outDir), progress/length(altFile)*100, "%"))
    }
    if (!is.null(ca)) {
        rownames(ca) <- partner_factors
        write.csv(ca, paste0(outDir, tf, '_test.csv'))
        gc(full=TRUE)
        rm(ca)
    }
}
factors <- list.files(intDir)
factors <- factors[factors %in% genes$name]
dir.create(oDir, showWarnings = FALSE, recursive = TRUE)
factors <- c("SALL2")
for (factor in factors) cobinding(factor, intDir, oDir)
#cl <- makeCluster(36, type = "PSOCK")
#registerDoParallel(cl)
#foreach(factor=factors, .packages = c("tidyverse", "data.table")) %dopar% cobinding(factor, intDir, oDir)
#stopCluster(cl)