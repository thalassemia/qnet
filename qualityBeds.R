library(tidyverse)

indexPath = "/home/sean/tfchip/human_factor_full_QC.txt"
dePath = "/home/sean/Research/SSvsP_SigTFs.csv"
bedPath = "/home/sean/tfchip/human_factor/"
outDir = "/home/sean/tfchip/goodBeds/"


index <- read_delim(indexPath, "\t") %>% filter(FastQC>=25 & UniquelyMappedRatio>=0.5  & PBC>=0.5 & PeaksUnionDHSRatio>=0.7)
# index <- read_delim(indexPath, "\t") %>% filter(Factor %in% c("DLX2", "DLX1", "DLX3", "ZNF695", "ATOH8", "MYB", 
#                                "NKX3-1", "MYBL2", "GL1", "CENPA", "SP6", 
#                                "ZNF850", "FOXM1", "ZNF492", "E2F1", "GATA2", 
#                                "E2F8", "E2F7", "KLF10", "FOSL1", "IRX6", 
#                                "PKNOX2", "EGR2", "IRF4", "NR3C2", "THRB", "OSR2",
#                                "TSHZ2", "KLF15", "MAF", "GRHL1", "BATF2", "EGR3",
#                                "SIX2", "ZNF540", "EBF4", "KLF4", "SATB1", "ZMAT1",
#                                "ZNF608"))
de <- read_csv(dePath) %>% rename(Factor = gene_name) %>% filter(padj <= 0.05)
factors <- inner_join(de, index, by="Factor") %>% select(DCid, Factor)

setwd(bedPath)
for (row in 1:nrow(factors)) {
  id = factors[row, "DCid"]
  tf = factors[row, "Factor"]
  for (file in list.files(pattern=paste(id, "_*", sep=""))) {
    dir.create(paste(outDir, tf, "/", sep=""), showWarnings=FALSE)
    file.copy(file, paste(outDir, tf, "/", id, ".bed", sep=""))
  }
}
