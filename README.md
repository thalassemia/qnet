# TFs and Quiescence
_**All code was designed to be run on the Hoffman2 cluster with 8 cores for Python scripts (~2GB RAM/core) and one core for R scripts (~4GB RAM) unless otherwise specified.**_

## Table of Contents
1. [Overview](#Overview)
2. [betaBatch.py](#betaBatchpy)
2. [cobinding.R](#cobindingR)
4. [deTargets.R](#deTargetsR)
5. [mergeDEandTF.R](#mergeDEandTFR)
6. [normalize.py](#normalizepy)
7. [peakOverlap.py](#peakOverlappy)
8. [qualityBeds.R](#qualityBedsR)
9. [sigTargets.py](#sigTargetspy)
10. [targetEnrichment.R](#targetEnrichmentR)


## [Overview](https://miro.com/app/board/o9J_km5-viQ=/)
![Workflow](workflow.png)


## [betaBatch.py](betaBatch.py)
**Dependencies:** Python 2.7 and [BETA](http://cistrome.org/BETA/)

Runs BETA minus (hg38) on input bed files to generate lists of potential targets ranked by regulatory potential as it is defined in [PMC4135175](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4135175/).


## [cobinding.R](cobinding.R)
_**Note:** Designed to run as a job array with 2 cores for each job (~4GB RAM/core)_

**Dependencies:** R 4.0+, tidyverse, pheatmap, dendsort, RColorBrewer, fastcluster

Creates dataframes combining all [rank normalized](#normalizepy), [peak overlapped](#peakOverlappy) bed files for each highly up/downregulated TF. After clustering using the incredibly memory efficient [fastcluster](http://danifold.net/fastcluster.html) library, it creates cobinding heatmaps as seen in [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/).


## [deTargets.R](deTargets.R)
**Dependencies:** R 4.0+ and tidyverse

Inner joins each [combined targets file](#sigTargetspy) with differential expression data generated using DESeq2.


## [mergeDEandTF.R](mergeDEandTF.R)
**Dependencies:** R 4.0+ and tidyverse

Performs an inner join between a [list of human transcription factors](http://humantfs.ccbr.utoronto.ca/download.php) and differential expression data generated using DESeq2. Sorts the resulting dataframe by ascending log2 fold change and removes any rows with:  

    abs(log2FoldChange) < 1.0 OR padj > 0.05


## [normalize.py](normalize.py)
**Dependencies:** Python 3.8+, tqdm, and pandas

Reads bed files, ranks peaks, and assigns them each a normalized rank score (see section C.2.2 in [supplementary paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/bin/NIHMS541492-supplement-Supplementary_Material.pdf) to [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/)) equal to `(R-r)/(R-1)`, where `r` is the rank of each peak in the total set of `R` peaks. Output files only retain columns for the chromosome, peak start, peak end, and normalized rank.

**Note:** Output files are ranked by descending signal value (enrichment) in `$SCRATCH/tfchip/normalize/signal/` and by ascending q value (p-adjusted) in `$SCRATCH/tfchip/normalize/qval/`.


## [peakOverlap.py](peakOverlap.py)
**Dependencies:** Python 3.8+, tqdm, pandas, and bedtools

For each transcription factor, runs:

    bedtools intersect -a {factor} -b {cofactor} -loj

`factor`: the focus TF's [highest quality](#qualityBedsR), [rank normalized](#normalizepy) bed file  
`cofactor`: any highest quality, rank normalized bed file (including the focus TF's)

Each row in the output files combines the chromosome, peak start, and peak end values of each peak in `factor` with the normalized rank of the overlapping peak in `cofactor` (-1 if no overlap; write multiple rows if more than one overlap).


## [qualityBeds.R](qualityBeds.R)
**Dependencies:** R 4.0+ and tidyverse

For each transcription factor that exhibits [significant differential expression](#mergeDEandTFR), this reads the [Cistrome ChIP-Seq](http://cistrome.org/db) index file to find relevant bed files. Then, it uses a mix of [ENCODE](https://www.encodeproject.org/data-standards/terms/) and [Cistrome](http://cistrome.org/db/#/about) ChIP-Seq quality guidelines to filter and rank all bed files. 

Specifically, bed files must at least meet the following criteria:

    UniquelyMappedRatio>=0.5 && fastQC>=25 && PBC>=0.5 && PeaksUnionDHSRatio>=0.7

The passing files are then sorted (first by descending `FRiP` then by descending `PeaksFoldChangeAbove10`) before being copied to a new directory with the following structure:

    $SCRATCH/goodBeds/{TF}/{rank}.bed

## [sigTargets.py](sigTargets.py)
**Dependencies:** Python 3.8+ and pandas

For each transcription factor, this consolidates all [BETA target predictions](#betaBatchpy) into a single tab-delimited csv file keeping the highest scoring instance of each duplicate target.


## [targetEnrichment.R](targetEnrichment.R)
**Dependencies:** R 4.0+, tidyverse, fastcluster, pheatmap, RColorBrewer

TFs are split into two groups, one for those that are upregulated with quiescence and one for those that are downregulated. All putative targets for all TFs in each group are consolidated at various [score](#betaBatchpy) thresholds (which represent BETA's confidence that a given target is a true target). 

_**Note:** Use as many cores as the desired number of thresholds to minimize run time (~4GB RAM/core)_

For each group and each score threshold, a matrix of `log2 Fold Changes` is compiled using the TFs in the group as rows and [all putative targets](#sigtargetspy) at that threshold as columns. A heatmap is generated using that matrix. Finally, a file is made listing ranks of the TFs when sorted by descending # of targets that go up with quiescence, # of targets that go down, ratio of # of up to # of downregulated targets, mean<sup>\*</sup> target enrichment with quiescence, 5<sup>th</sup> %ile<sup>\*</sup> target enrichment, and 95<sup>th</sup> %ile<sup>\*</sup> target enrichment.

Outputs in a new directory with the following naming scheme:

    $SCRATCH/output/enrich/{updown}{threshold}{type}

`updown`: "Up" for group of TFs upregulated with quiescence, "Down" for group of TFs downregulated with quiescence  
`threshold`: score threshold (see regulatory potential formula in [PMC4135175](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4135175/) for more details) used to subset targets  
`type`: "count" for target count statistics, "matrix" for enrichment matrices, or "heatmap" for heatmaps

<sup>\*</sup> Zeros are excluded in these calculations as they represent non-targets for any given TF