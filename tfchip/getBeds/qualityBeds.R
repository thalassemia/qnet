library(tidyverse)

# input and output files and directories
indexPath = system("echo $SCRATCH/tfchip/human_factor_full_QC.txt", intern = TRUE)
deTFPath = system("echo $SCRATCH/data/deTF.csv", intern = TRUE)
bedPath = system("echo $SCRATCH/tfchip/human_factor/", intern = TRUE)
outDir = system("echo $SCRATCH/cobinding/beds/", intern = TRUE)
# tfs for which there is already data
existingFiles <- trimws(list.files(outDir))
encTFs = data.frame(DCid = rep('N/A', length(existingFiles)), Factor = existingFiles, rank = rep('N/A', length(existingFiles)), 
                    Cell_line = rep('N/A', length(existingFiles)), Cell_type = rep('N/A', length(existingFiles)), Tissue_type = rep('N/A', length(existingFiles)),
                    FastQC = rep('N/A', length(existingFiles)), UniquelyMappedRatio = rep('N/A', length(existingFiles)),
                    PBC = rep('N/A', length(existingFiles)), PeaksFoldChangeAbove10 = rep('N/A', length(existingFiles)),
                    FRiP = rep('N/A', length(existingFiles)), PeaksUnionDHSRatio = rep('N/A', length(existingFiles)))

# find files matching minimum quality criteria and assign them a rank by descending FRiP
index <- read_delim(indexPath, "\t") %>% filter(FastQC>=25 & UniquelyMappedRatio>=0.5  & PBC>=0.5)
index <- arrange(index, desc(FRiP), desc(PeaksFoldChangeAbove10), desc(PeaksUnionDHSRatio))
index$rank <- rownames(index)

# figure out which TFs have corresponding bed files
de <- read_csv(deTFPath) %>% rename(Factor = gene_name)
cistromeTFs <- inner_join(index, de, by="Factor") %>% select(DCid, Factor, rank, Cell_line, Cell_type, Tissue_type, FastQC, UniquelyMappedRatio, PBC, PeaksFoldChangeAbove10, FRiP, PeaksUnionDHSRatio, log2FoldChange, padj)
encodeTFs <- inner_join(encTFs, de, by="Factor") %>% select(DCid, Factor, rank, Cell_line, Cell_type, Tissue_type, FastQC, UniquelyMappedRatio, PBC, PeaksFoldChangeAbove10, FRiP, PeaksUnionDHSRatio, log2FoldChange, padj)
rbind(cistromeTFs, encodeTFs) %>% write_csv(paste(outDir, "key.csv", sep=""))
anti_join(de, index, by="Factor") %>% anti_join(encTFs, by = "Factor") %>% select(Factor) %>% write_csv(paste(outDir, "nodata.csv", sep=""))

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
