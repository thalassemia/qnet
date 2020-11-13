# Create Co-Binding Heatmaps

## Table of Contents
1. [TF-centric Approach](#tf-centric-approach)
1. [DE Gene Approach](#de-gene-approach)

## [TF-centric Approach](cobindO.R)
_**Note:** Designed to run as a job array with 2 cores for each job (~4GB RAM/core)_

**Dependencies:** R 4.0+, tidyverse, pheatmap, dendsort, RColorBrewer, fastcluster

Creates separate dataframes for each highly enriched or depleted TF by concatenating all the rank-normalized, peak-overlapped bed files for that individual TF. After performing hierarchical clustering on both the rows and columns (Euclidean distance, Ward linkage) using the incredibly memory efficient [fastcluster](http://danifold.net/fastcluster.html) library, it creates cobinding heatmaps as seen in [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/). 

For each heatmap, the titular TF is the so-called "focus factor" for that graphic. The value at a given row (potential cobinding factor) and column (focus factor peak) is the rank normalized score of the potential cobinding factor peak that overlaps with the focus factor peak. If no overlapping peak exists, a score of 0 is assigned. If multiple overlapping peaks exist, the highest score among them is assigned.

By default, heatmaps are created for each TF in the top 20 by log2 fold change DE for which there is data. The rows of the heatmaps are all DE TF's that change in the same direction as the focus factor (all enriched or all depleted with quiescence).

## [DE Gene Approach](cobindU.r)

**Dependencies:** R 4.0+, tidyverse, pheatmap, dendsort, RColorBrewer, fastcluster

Creates a heatmap visualization of an already consolidated co-binding matrix (see [Preprocessing Data](../preprocess/info.md)).