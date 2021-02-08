# Histone Modifications

## Table of Contents
1. [Filter ChIP-seq Data](#filter-chip-seq-data)
1. [Annotate Peaks](#annotate-peaks)
1. [Preprocess Data](#preprocess-data)
1. [Create Co-binding Map](#create-co-binding-maps)

---

## [Filter ChIP-seq Data](qualityBeds.R)
Refer to [Cistrome ChIP-seq Data](../tfchip/info.md#Cistrome-chip-seq-data) for more details. The only notable change here is that all histone bed files (as opposed to only DE TF bed files in the link above) passing the quality criteria are ranked and copied into the new directory structure.

## [Annotate Peaks](annotate.R)
For each bed file in a given directory, assigns each peak to its nearest transcript annotation using ChIPseeker (promoter set to TSS ± 3000bp) and GENCODE V29.

## [Preprocess Data](preprocess.py)
Refer to [DE Gene Prep](../tfchip/info.md#de-gene-prep) for more details. The only notable change here is that in step 2, one can choose to either filter for peaks annotated to promoters of significant DE genes (RNA-seq, DEseq2) OR significant DA regions (ATAC-seq, DiffBind). In either case, significance refers to having an adjusted p-value less than 0.05.

## [Create Co-binding Map](cobinding.R)
Creates a heatmap visualization of all co-binding matrices in a given directory. For each heatmap, hierarchical clustering is performed on both the rows and columns using Euclidean distance and Ward linkage.