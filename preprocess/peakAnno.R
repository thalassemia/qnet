library(doParallel)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

dir <- system("echo $SCRATCH/tfchip/lisabeds/", intern = TRUE)
setwd(dir)
baseOutDir <- system("echo $SCRATCH/tfchip/lisaannotated/", intern = TRUE)

factors <- list.files()
factors <- factors[!(factors %in% c("nodata.csv", "key.csv"))] 

print("Starting annotation")

cl <- makeCluster(36, type = "PSOCK")
registerDoParallel(cl)
factors <- foreach(factor = factors, .packages = c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene")) %dopar% {
  setwd(file.path(dir, factor))
  beds <- list.files()
  outDir <- file.path(baseOutDir, factor)
  dir.create(outDir, showWarnings = FALSE, recursive = T)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  
  for(bed in beds){  
    setwd(file.path(dir, factor))
    peak <- ChIPseeker::readPeakFile(bed)
    suppressWarnings(peakAnno <- ChIPseeker::annotatePeak(peak, TxDb = txdb, annoDb = "org.Hs.eg.db"))
    peakAnno <- as.data.frame(peakAnno)
    setwd(outDir)
    write.csv(peakAnno, paste0(bed,".csv"), row.names = F)
  } 
  factor
} 
stopCluster(cl)