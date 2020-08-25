# TFs and Quiescence


## Table of Contents
1. [Overview](#Overview)
2. [betaBatch.py](#betaBatchpy)
2. [cobinding.R](#cobindingR)
4. [deTargets.R](#deTargetsR)
5. [mergeDEandTF.R](#mergeDEandTFR)
6. [normalize.py](#normalizepy)
7. [peakOverlap.py](#peakOverlappy)
8. [pullPeakFiles.py](#pullPeakFilespy)
9. [sigTargets.py](#sigTargetspy)
10. [targetEnrichment.R](#targetEnrichmentR)


## [Overview](https://miro.com/app/board/o9J_kmhzx0k=/?moveToWidget=3074457349492140638&cot=12)
![Workflow](workflow.png)


## [betaBatch.py](betaBatch.py)
**Dependencies:** Python 2.7 and [BETA](http://cistrome.org/BETA/)

Runs BETA minus (hg38) on input bed files to generate lists of potential targets ranked by regulatory potential as it is defined in [PMC4135175](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4135175/).


## [cobinding.R](cobinding.R)
_**Note:** Designed to run as a job array on Hoffman2 Cluster_

**Dependencies:** R 4.0+, tideverse, pheatmap, dendsort, RColorBrewer, fastcluster

Creates dataframes combining all [rank normalized](#normalizepy), [peak overlapped](#peakOverlappy) bed files for each highly up/downregulated TF. After clustering using the incredibly memory efficient [fastcluster](http://danifold.net/fastcluster.html) library, it creates cobinding heatmaps as seen in [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/).


## [deTargets.R](deTargets.R)
**Dependencies:** R 4.0+ and tidyverse

Inner join each [combined targets file](#sigTargetspy) with the differential expression data.


## [mergeDEandTF.R](mergeDEandTF.R)
**Dependencies:** R 4.0+ and tidyverse

Performs an inner join between a [list of human transcription factors](http://humantfs.ccbr.utoronto.ca/download.php) and differential expression data generated using DESeq2. Sorts the resulting dataframe by ascending log2 fold change and removes any rows with:  
>`abs(log2FoldChange)<1.0 || padj>0.05`


## [normalize.py](normalize.py)
**Dependencies:** Python 3.8+, tqdm, and pandas

Reads bed files, ranks<sup>*</sup> peaks, and assigns them each a normalized rank score (see section C.2.2 in [supplementary paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/bin/NIHMS541492-supplement-Supplementary_Material.pdf) to [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/)) equal to `(R-r)/(R-1)`, where `r` is the rank of each peak in the total set of `R` peaks. Output files only retain columns for the chromosome, peak start, peak end, and normalized rank.

<sup>*</sup> Peaks are ranked by descending signal value (enrichment) when `rankBy == 6` and by ascending q value (p-adjusted) when `rankBy == 8`.


## [peakOverlap.py](peakOverlap.py)
**Dependencies:** Python 3.8+, tqdm, pandas, and bedtools

For each transcription factor, runs all unique combinations of:
>`bedtools intersect -a {factor} -b {cofactor} -loj` 

`factor` and `cofactor` are each a given TF's [highest quality](#qualityBedsR), [rank normalized](#normalizepy) bed file.  

For each unique combination of `factor` and `cofactor`, each row in the output file combines the chromosome, peak start, and peak end values from `factor` with the normalized rank of the overlapping peaks in `cofactor` (-1 if none; multiple rows if more than one).


## [pullPeakFiles.py](pullPeakFiles.py)
**Dependencies:** Python 3.8+ and pandas

For each transcription factor that exhibits [significant differential expression](#mergeDEandTFR), reads [Cistrome ChIP-Seq](http://cistrome.org/db) index file to copy TF-associated narrowPeak bed files to a new directory with the following structure:
>`~/peaks/{factor}/*.bed`


## [qualityBeds.R](qualityBeds.R)
**Dependencies:** R 4.0+ and tidyverse

Uses a mix of [ENCODE](https://www.encodeproject.org/data-standards/terms/) and [Cistrome](http://cistrome.org/db/#/about) ChIP-Seq quality guidelines to filter and rank the quality of all bed files in a new directory with the following structure:
>`~/goodBeds/{factor}/{rank}.bed`

Specifically, bed files have to meet the following criteria:
>`UniquelyMappedRatio>=0.5 && fastQC>=25 && PBC>=0.5 && PeaksUnionDHSRatio>=0.7`

The passing files were then sorted first by descending `FRiP` then by descending `PeaksFoldChangeAbove10` before being assigned their final ranks.


## [sigTargets.py](sigTargets.py)
**Dependencies:** Python 3.8+ and pandas

For each transcription factor, consolidates all BETA output files into a single tab-delimited csv file with all lower scoring duplicate targets dropped.


## [targetEnrichment.R](targetEnrichment.R)
_**Note:** Very memory intensive, 32G on Hoffman2_

**Dependencies:** R 4.0+, tidyverse, fastcluster, pheatmap, RColorBrewer

Create matrix of `log2 Fold Change * regulatory potential` values with all putative targets as columns and TFs as rows. Plot heatmap. Create lists of TFs by most putative targets up/down with quiescence and highest/lowest mean value.