library(tidyverse)

indexPath = system("echo $SCRATCH/hmchip/human_hm_full_QC.txt", intern = TRUE)
bedPath = system("echo $SCRATCH/hmchip/human_hm/", intern = TRUE)
outDir = system("echo $SCRATCH/hmchip/sorted/", intern = TRUE)
dir.create(outDir)

index <- read_delim(indexPath, "\t") %>% filter(FastQC>=25 & UniquelyMappedRatio>=0.5  & PBC>=0.5)
index <- arrange(index, desc(FRiP), desc(PeaksFoldChangeAbove10), desc(PeaksUnionDHSRatio))
index$rank <- rownames(index)

for (row in 1:nrow(index)) {
  id = index[row, "DCid"]
  rank = index[row, "rank"]
  hm = index[row, "Factor"]
  for (file in list.files(bedPath, pattern=paste(id, "_*", sep=""))) {
    dir.create(paste(outDir, hm, "/", sep=""), showWarnings=FALSE)
    file.copy(paste(bedPath, file, sep=""), paste(outDir, hm, "/", rank, ".bed", sep=""), overwrite=TRUE)
  }
}