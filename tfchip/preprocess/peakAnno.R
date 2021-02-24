library(doParallel)
library(tidyverse)
library(ChIPseeker)
library(org.Hs.eg.db)

dir <- system("echo $SCRATCH/cobinding/beds/", intern = TRUE)
setwd(dir)
baseOutDir <- system("echo $SCRATCH/cobinding/annotated/", intern = TRUE)

factors <- list.files()
# exclude informational files
factors <- factors[!(factors %in% c("encodeKey.csv", "nodata.csv", "key.csv"))] 

print("Starting annotation")

# 36 threads, ~70 GB memory
cl <- makeCluster(36, type = "PSOCK")
registerDoParallel(cl)
factors <- foreach(factor = factors, .packages = c("ChIPseeker", "org.Hs.eg.db")) %dopar% {
  setwd(file.path(dir, factor))
  beds <- list.files()
  outDir <- file.path(baseOutDir, factor)
  dir.create(outDir, showWarnings = FALSE, recursive = T)

  # load GENCODE V29 database
  annoFile <- system("echo $SCRATCH/data/gencode_v29.sqlite", intern = TRUE)
  txdb <- loadDb(annoFile)
  
  # annotate and output to new directory
  for(bed in beds){  
    setwd(file.path(dir, factor))
    peak <- ChIPseeker::readPeakFile(bed)
    suppressWarnings(peakAnno <- ChIPseeker::annotatePeak(peak, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE))
    peakAnno <- as.data.frame(peakAnno)
    setwd(outDir)
    write.csv(peakAnno, bed, row.names = F)
  } 
  factor
} 
stopCluster(cl)