# Filter ChIP-Seq Data

## Table of Contents
1. [Get DE TF List](#get-de-tf-list)
1. [ChEA3/Lisa Output](#chea3/lisa-output-chip-seq-data)
1. [DE TFs ChIP-seq Data](#de-tfs-chip-seq-data)

## [Get DE TF List](mergeDEandTF.R)
**Dependencies:** R 4.0+ and tidyverse

Performs an inner join between a [list of human transcription factors](http://humantfs.ccbr.utoronto.ca/download.php) and differential expression data generated using DESeq2. Sorts the resulting dataframe by ascending log2 fold change and removes any rows with:  

    abs(log2FoldChange) < 1.0 OR padj > 0.05

## [ChEA3/Lisa Output ChIP-Seq Data](chea_lisa_beds.r)
**Dependencies:** R 4.0+ and tidyverse

The top 500 DE genes were input into [ChEA3](https://maayanlab.cloud/chea3/) and [Lisa](http://lisa.cistrome.org/) to generate two distinct lists of likely transcriptional regulators of quiescence. For each of the top 200 TFs from each list, this reads the the [Cistrome ChIP-Seq](http://cistrome.org/db) index file to find relevant bed files. Then, it uses a mix of [ENCODE](https://www.encodeproject.org/data-standards/terms/) and [Cistrome](http://cistrome.org/db/#/about) ChIP-Seq quality guidelines to filter and rank all bed files. 

Specifically, bed files must at least meet the following criteria:

    UniquelyMappedRatio>=0.5 && fastQC>=25 && PBC>=0.5 && PeaksUnionDHSRatio>=0.7

The passing files are then sorted (first by descending `FRiP` then by descending `PeaksFoldChangeAbove10`) before being copied to a new directory with the following structure:

    $SCRATCH/{chea OR lisa}Beds/{TF}/{rank}.bed

## [DE TFs ChIP-Seq Data](qualityBeds.R)
**Dependencies:** R 4.0+ and tidyverse

For each transcription factor that exhibits [significant differential expression](#mergeDEandTFR), this reads the [Cistrome ChIP-Seq](http://cistrome.org/db) index file to find relevant bed files. Then, it uses a mix of [ENCODE](https://www.encodeproject.org/data-standards/terms/) and [Cistrome](http://cistrome.org/db/#/about) ChIP-Seq quality guidelines to filter and rank all bed files. 

The passing files are then sorted (first by descending `FRiP` then by descending `PeaksFoldChangeAbove10`) before being copied to a new directory with the following structure:

    $SCRATCH/goodBeds/{TF}/{rank}.bed