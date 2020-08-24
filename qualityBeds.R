library(tidyverse)

indexPath = "/u/scratch/s/seanchea/human_factor_full_QC.txt"
dePath = "/u/scratch/s/seanchea/deTF.csv"
bedPath = "/u/scratch/s/seanchea/tfbed/human_factor/"
outDir = "/u/scratch/s/seanchea/goodBeds/"
dir.create(outDir)

index <- read_delim(indexPath, "\t") %>% filter(FastQC>=25 & UniquelyMappedRatio>=0.5  & PBC>=0.5 & PeaksUnionDHSRatio>=0.7)
index <- arrange(index, desc(FRiP), desc(PeaksFoldChangeAbove10))
index$rank <- rownames(index)
# index <- read_delim(indexPath, "\t") %>% filter(Factor %in% c("DLX2", "DLX1", "DLX3", "ZNF695", "ATOH8", "MYB", 
#                                "NKX3-1", "MYBL2", "GL1", "CENPA", "SP6", 
#                                "ZNF850", "FOXM1", "ZNF492", "E2F1", "GATA2", 
#                                "E2F8", "E2F7", "KLF10", "FOSL1", "IRX6", 
#                                "PKNOX2", "EGR2", "IRF4", "NR3C2", "THRB", "OSR2",
#                                "TSHZ2", "KLF15", "MAF", "GRHL1", "BATF2", "EGR3",
#                                "SIX2", "ZNF540", "EBF4", "KLF4", "SATB1", "ZMAT1",
#                                "ZNF608"))
de <- read_csv(dePath) %>% rename(Factor = gene_name) %>% filter(padj <= 0.05)
factors <- inner_join(index, de, by="Factor") %>% select(DCid, Factor, rank)
factors %>% write_csv(paste(outDir, "key.csv", sep=""), append=FALSE)

for (row in 1:nrow(factors)) {
  id = factors[row, "DCid"]
  rank = factors[row, "rank"]
  tf = factors[row, "Factor"]
  for (file in list.files(bedPath, pattern=paste(id, "_s*", sep=""))) {
    dir.create(paste(outDir, tf, "/", sep=""), showWarnings=FALSE)
    file.copy(paste(bedPath, file, sep=""), paste(outDir, tf, "/", rank, ".bed", sep=""), overwrite=TRUE)
  }
}
