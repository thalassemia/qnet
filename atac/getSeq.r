library(tidyverse)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

Hsapiens <- BSgenome.Hsapiens.UCSC.hg38

atac <- read_csv("/u/scratch/s/seanchea/atac/atacDE.csv")
ranges <- GRanges(
    seqnames = Rle(atac$seqnames),
    ranges = IRanges(atac$start, end = atac$end)
)

seq <- getSeq(Hsapiens,ranges)
names(seq) <- seq_len(length(ranges))
writeXStringSet(seq,file="/u/scratch/s/seanchea/atac/DEsequences")