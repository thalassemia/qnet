library(tidyverse)

# input and output files and directories
indexPath = system("echo $SCRATCH/tfchip/human_factor_full_QC.txt", intern = TRUE)
deTFPath = system("echo $SCRATCH/data/deTF.csv", intern = TRUE)
bedPath = system("echo $SCRATCH/tfchip/human_factor/", intern = TRUE)
outDir = system("echo $SCRATCH/cobinding/beds/", intern = TRUE)
# tfs for which there is already data
existingFiles <- trimws(list.files(outDir))
enc <- read_csv(paste0(outDir, "encodeKey.csv")) %>% arrange(desc(FRiP))
encTFs <- data.frame(DCid = enc[['X1']], Factor = enc[['Target']], rank = rownames(enc), 
                    Cell_line = enc[['Cell Type(s)']], Cell_type = rep('N/A', nrow(enc)), Tissue_type = rep('N/A', nrow(enc)),
                    FastQC = rep('N/A', nrow(enc)), UniquelyMappedRatio = enc[['NRF']],
                    PBC = enc[['PBC1']], PeaksFoldChangeAbove10 = rep('N/A', nrow(enc)),
                    FRiP = enc[['FRiP']], PeaksUnionDHSRatio = rep('N/A', nrow(enc)))

# find files matching minimum quality criteria and assign them a rank by descending FRiP
index <- read_delim(indexPath, "\t") %>% filter(FastQC>=25 & UniquelyMappedRatio>=0.6  & PBC>=0.8 & PeaksUnionDHSRatio>=0.7)
index <- arrange(index, desc(FRiP), desc(PeaksFoldChangeAbove10))
index$rank <- rownames(index)

# figure out which TFs have corresponding bed files
de <- read_csv(deTFPath) %>% rename(Factor = gene_name)
cistromeTFs <- inner_join(index, de, by="Factor") %>% select(DCid, Factor, rank, Cell_line, Cell_type, Tissue_type, FastQC, UniquelyMappedRatio, PBC, PeaksFoldChangeAbove10, FRiP, PeaksUnionDHSRatio, log2FoldChange, padj)
encodeTFs <- inner_join(encTFs, de, by="Factor") %>% select(DCid, Factor, rank, Cell_line, Cell_type, Tissue_type, FastQC, UniquelyMappedRatio, PBC, PeaksFoldChangeAbove10, FRiP, PeaksUnionDHSRatio, log2FoldChange, padj)
rbind(cistromeTFs, encodeTFs) %>% write_csv(paste(outDir, "key.csv", sep=""))
anti_join(de, index, by="Factor") %>% anti_join(encTFs, by = "Factor") %>% select(Factor) %>% write_csv(paste(outDir, "nodata.csv", sep=""))
# get count of downregulated DE TFs with data
print(paste("Down:", unique(filter(rbind(cistromeTFs, encodeTFs), log2FoldChange < 0)[, 'Factor'])))
# get count of upregulated DE TFs with data
print(paste("Up:", unique(filter(rbind(cistromeTFs, encodeTFs), log2FoldChange >= 0)[, 'Factor'])))

for (row in 1:nrow(cistromeTFs)) {
  id = cistromeTFs[row, "DCid"]
  rank = cistromeTFs[row, "rank"]
  tf = cistromeTFs[row, "Factor"]
  # copy all relevant files to a folder named the corresponding TF 
  for (file in list.files(bedPath, pattern=paste(id, "_sort*", sep=""))) {
    dir.create(paste(outDir, tf, "/", sep=""), showWarnings=FALSE)
    file.copy(paste(bedPath, file, sep=""), paste(outDir, tf, "/", rank, ".bed", sep=""), overwrite=TRUE)
  }
}
