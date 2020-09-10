library(tidyverse)
library(rGREAT)
library(GenomicRanges)

peak_dir <- system("echo $SCRATCH/goodBeds", intern = T)
out_dir <- system("echo $SCRATCH/GO", intern = T)

factors <- list.files(peak_dir)
factors <- factors[!factors %in% c("nodata.csv", "key.csv")]

for (factor in factors) {
    setwd(file.path(peak_dir, factor))
    beds <- list.files()
    dir.create(file.path(out_dir, factor), showWarnings = F, recursive = T)
    for (bed in beds) {
        path <- file.path(peak_dir, factor, bed)
        df <- read.csv(path, header = F, sep = "\t")
        df <- GRanges(seqnames = df$V1,
                      ranges = IRanges(start = df$V2,
                                       end = df$V3,
                                       names = df$V4))
        new_style <- seqlevels(df) %>% mapSeqlevels("UCSC")
        new_style <- new_style[complete.cases(new_style)]
        df <- renameSeqlevels(df, new_style)
        werid_names <- seqlevels(df)[c(grep("chrGL*", seqlevels(df)),
                        grep("chrKI*", seqlevels(df)),
                        grep("chrMT*", seqlevels(df)))]
        df <- dropSeqlevels(df, werid_names, pruning.mode = "coarse")
        job <- submitGreatJob(df, species = "hg38")
        print(job)
        data <- getEnrichmentTables(job, download_by = "tsv")
        res <- plotRegionGeneAssociationGraphs(job)
        out_path <- file.path(out_dir, factor, bed)
        write.csv(data[[1]], paste0(out_path, "_GOBP.csv"), row.names = F)
        write.csv(data[[2]], paste0(out_path, "_GOCC.csv"), row.names = F)
        write.csv(data[[3]], paste0(out_path, "_GOMF.csv"), row.names = F)
        write.csv(res, paste0(out_path, "_genes.csv"), row.names = F)
        print(paste0("Done with ", bed))
    }
}
