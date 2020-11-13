library(tidyverse)

indexPath = system("echo $SCRATCH/tfchip/human_factor_full_QC.txt", intern = TRUE)
cheaPath = system("echo $SCRATCH/data/chea3.tsv", intern = TRUE)
bedPath = system("echo $SCRATCH/tfchip/human_factor/", intern = TRUE)
outDir = system("echo $SCRATCH/tfchip/cheabeds/", intern = TRUE)
dir.create(outDir)

chea <- read_delim(cheaPath, "\t")
chea <- rename(chea, Factor = TF)
chea <- chea[1:200,]

index <- read_delim(indexPath, "\t") %>% filter(FastQC>=25 & UniquelyMappedRatio>=0.5  & PBC>=0.5)
index <- arrange(index, desc(FRiP), desc(PeaksFoldChangeAbove10), desc(PeaksUnionDHSRatio))
index$rank <- rownames(index)

factors <- inner_join(index, chea, by="Factor") %>% select(DCid, Factor, rank, Cell_line, Cell_type, Tissue_type, FastQC, UniquelyMappedRatio, PBC, PeaksFoldChangeAbove10, FRiP, PeaksUnionDHSRatio)
anti_join(chea, index, by="Factor") %>% select(Factor) %>% write_csv(paste(outDir, "nodata.csv", sep=""))
factors %>% write_csv(paste(outDir, "key.csv", sep=""))

for (row in 1:nrow(factors)) {
  id = factors[row, "DCid"]
  rank = factors[row, "rank"]
  tf = factors[row, "Factor"]
  for (file in list.files(bedPath, pattern=paste(id, "_s*", sep=""))) {
    dir.create(paste(outDir, tf, "/", sep=""), showWarnings=FALSE)
    file.copy(paste(bedPath, file, sep=""), paste(outDir, tf, "/", rank, ".bed", sep=""), overwrite=TRUE)
  }
}

lisaPath = system("echo $SCRATCH/data/lisa.csv", intern = TRUE)
outDir = system("echo $SCRATCH/tfchip/lisabeds/", intern = TRUE)
dir.create(outDir)

lisa <- read_csv(lisaPath) %>% select("Transcription Factor")
lisa <- rename(lisa, Factor = "Transcription Factor")
lisa <- lisa[1:200,]

factors <- inner_join(index, lisa, by="Factor") %>% select(DCid, Factor, rank, Cell_line, Cell_type, Tissue_type, FastQC, UniquelyMappedRatio, PBC, PeaksFoldChangeAbove10, FRiP, PeaksUnionDHSRatio)
anti_join(lisa, index, by="Factor") %>% select(Factor) %>% write_csv(paste(outDir, "nodata.csv", sep=""))
factors %>% write_csv(paste(outDir, "key.csv", sep=""))

for (row in 1:nrow(factors)) {
  id = factors[row, "DCid"]
  rank = factors[row, "rank"]
  tf = factors[row, "Factor"]
  for (file in list.files(bedPath, pattern=paste(id, "_s*", sep=""))) {
    dir.create(paste(outDir, tf, "/", sep=""), showWarnings=FALSE)
    file.copy(paste(bedPath, file, sep=""), paste(outDir, tf, "/", rank, ".bed", sep=""), overwrite=TRUE)
  }
}