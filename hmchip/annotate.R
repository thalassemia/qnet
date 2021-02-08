library(doParallel)

dir <- system("echo $SCRATCH/hmchip/sorted/", intern = TRUE)
setwd(dir)
baseOutDir <- system("echo $SCRATCH/hmchip/annotated/", intern = TRUE)

factors <- list.files()
factors <- factors[!(factors %in% c("nodata.csv", "key.csv"))] 

print("Starting annotation")

cl <- makeCluster(36, type = "PSOCK")
registerDoParallel(cl)
factors <- foreach(factor = factors, .packages = c("ChIPseeker", "GenomicFeatures")) %dopar% {
  # load GENCODE V29 database
  annoFile <- system("echo $SCRATCH/data/gencode_v29.sqlite", intern = TRUE)
  txdb <- loadDb(annoFile)
  
  setwd(file.path(dir, factor))
  beds <- list.files()
  # only annotate the top-ranked bed file for each mark to save time
  beds <- beds[order(as.integer(sapply(beds, function(x) strsplit(x, ".")[[1]][1], simplify=TRUE)), decreasing=FALSE)]
  bed <- beds[1]
  outDir <- file.path(baseOutDir, factor)
  dir.create(outDir, showWarnings = FALSE, recursive = T)
  peak <- ChIPseeker::readPeakFile(bed)
  suppressWarnings(peakAnno <- ChIPseeker::annotatePeak(peak, TxDb = txdb, annoDb = "org.Hs.eg.db"))
  peakAnno <- as.data.frame(peakAnno)
  setwd(outDir)
  write.csv(peakAnno, paste0(bed,".csv"), row.names = F) 
  factor
} 
stopCluster(cl)