library(tidyverse)

indexPath = system("echo $SCRATCH/tfchip/human_factor_full_QC.txt", intern = TRUE)
deTFPath = system("echo $SCRATCH/data/deTF.csv", intern = TRUE)
bedPath = system("echo $SCRATCH/tfchip/human_factor/", intern = TRUE)
outDir = system("echo $SCRATCH/tfchip/goodBeds/", intern = TRUE)
dir.create(outDir)

index <- read_delim(indexPath, "\t") %>% filter(FastQC>=25 & UniquelyMappedRatio>=0.5  & PBC>=0.5 & PeaksUnionDHSRatio>=0.7)
index <- arrange(index, desc(FRiP), desc(PeaksFoldChangeAbove10))
index$rank <- rownames(index)

de <- read_csv(deTFPath) %>% rename(Factor = gene_name) %>% filter(padj <= 0.05)
factors <- inner_join(index, de, by="Factor") %>% select(DCid, Factor, rank)
anti_join(de, index, by="Factor") %>% select(Factor) %>% write_csv(paste(outDir, "nodata.csv", sep=""))
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
