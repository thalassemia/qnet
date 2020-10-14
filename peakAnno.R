library(doParallel)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

dir <- system("echo $SCRATCH/tfchip/goodBedsSALL2/", intern = TRUE)
setwd(dir)
baseOutDir <- system("echo $SCRATCH/tfchip/goodAnnoSALL2/", intern = TRUE)

factors <- list.files()
factors <- factors[!(factors %in% c("nodata.csv", "key.csv"))] 

print("Starting annotation")

cl <- makeCluster(30, type = "PSOCK")
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
    suppressWarnings(peakAnno <- ChIPseeker::annotatePeak(peak, TxDb = txdb))
    peakAnno <- as.data.frame(peakAnno)
    setwd(outDir)
    write.csv(peakAnno, paste0(bed,".csv"), row.names = F)
  } 
  factor
}  