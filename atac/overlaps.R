library(GenomicRanges)
library(regioneR)
library(data.table)

DNasePath <- system("echo $SCRATCH/data/roadmapDNase10.bed", intern = TRUE)
#DNaseI <- toGRanges(DNasePath)
#print(head(DNaseI))
atacPath <- system("echo $SCRATCH/data/SS_P_diff_bound_sites.csv", intern = TRUE)
atac <- fread(atacPath)
#print(head(atac))
outPath <- system("echo $SCRATCH/atac/merged_with_roadmap_DHS.bed", intern = TRUE)

tempATAC <- tempfile()
write.table(atac,file=tempATAC,quote=F,sep="\t",col.names=F,row.names=F)

system(paste("bedtools intersect -wao -a", tempATAC, "-b", DNasePath, ">", outPath))

chromHMMPath <- system("echo $SCRATCH/data/E126_15_coreMarks_hg38lift_mnemonics.bed", intern = TRUE)
#chromHMM <- toGRanges(chromHMMPath)
#print(head(chromHMM))
outPath <- system("echo $SCRATCH/atac/merged_with_roadmap_chromHMM.bed", intern = TRUE)
system(paste("bedtools intersect -wao -a", tempATAC, "-b", chromHMMPath, ">", outPath))

ccrePath <- system("echo $SCRATCH/data/GRCh38-ccREs.bed", intern = TRUE)
outPath <- system("echo $SCRATCH/atac/merged_with_encode_ccre.bed", intern = TRUE)
system(paste("bedtools intersect -wao -a", tempATAC, "-b", ccrePath, ">", outPath))
unlink(tempATAC)