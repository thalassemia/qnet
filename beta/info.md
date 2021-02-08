# Get DE TF Targets

## Table of Contents
1. [betaBatch.py](#betaBatchpy)
1. [sigTargets.py](#sigTargetspy)
1. [deTargets.R](#deTargetsR)
1. [targetEnrichment.R](#targetEnrichmentR)

## [betaBatch.py](betaBatch.py)
**Dependencies:** Python 2.7 and [BETA](http://cistrome.org/BETA/)

Runs BETA minus (hg38) on input bed files to generate lists of potential targets ranked by regulatory potential as it is defined in [PMC4135175](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4135175/).

## [sigTargets.py](sigTargets.py)
For each transcription factor, this consolidates all [BETA target predictions](#betaBatchpy) into a single tab-delimited csv file keeping the highest scoring instance of each duplicate target.

## [deTargets.R](deTargets.R)
Inner joins each [combined targets file](#sigTargetspy) with differential expression data generated using DESeq2.

## [targetEnrichment.R](targetEnrichment.R)
_**Note:** Use 7 cores (~4GB RAM/core)_

TFs are split into two groups, one for those that are upregulated with quiescence and another for those that are downregulated. For each group, a list of unique putative targets (meaning those that show up in the BETA output of at least one TF in the group) is compiled at score thresholds of 0, 0.5, 1.0, 1.5, 2.0, 2.5, and 3.0. 

For each group and score threshold, this creates:
* A "target enrichment" matrix. The value displayed at any given row (TF) and column (putative target) is the SSvsP log2 fold change (l2FC) of the putative target, defaulting to 0 if the target was not predicted for the TF.
* A heatmap of the above matrix after performing hierarchical clustering on both the rows and columns using Euclidean distance and Ward linkage. Note that the color gradient is not linear but rather centered around zero and partitioned by percentile.
* A file containing several statistics of interest. For each TF within the group, this tabulates the # of targets that go up with quiescence, # of targets that go down with quiescence, ratio of # of targets that go up to # that go down, mean l2FC, 5<sup>th</sup> %ile l2FC, and 95<sup>th</sup> %ile l2FC. Note that all zeros were excluded when determing means and percentiles. This file also shows the ranks associated with each TF when sorting (descending) by each of the aforementioned statistics. 

Outputs in a new directory with the following naming scheme:

    $SCRATCH/output/enrich/{updown}{threshold}{type}

`updown`: "Up" for group of TFs that are upregulated with quiescence, "Down" for group of TFs that are downregulated with quiescence  
`threshold`: score threshold (see regulatory potential formula in [PMC4135175](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4135175/) for more details) used to subset targets  
`type`: "matrix" for target enrichment matrices, "heatmap" for heatmaps, or "count" for statistics of interest