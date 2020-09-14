# TFs and Quiescence
_**All code was designed to be run on the Hoffman2 cluster with 8 cores for Python scripts (~2GB RAM/core) and one core for R scripts (~4GB RAM) unless otherwise specified.**_

## Table of Contents
1. [Overview](#Overview)
1. [betaBatch.py](#betaBatchpy)
1. [cobinding.R](#cobindingR)
1. [deTargets.R](#deTargetsR)
1. [fibroblast.py](#fibroblastpy)
1. [greatBatch.R](#greatbatchR)
1. [H4K20me3.py](#H4K20me3py)
1. [mergeDEandTF.R](#mergeDEandTFR)
1. [normalize.py](#normalizepy)
1. [peakOverlap.py](#peakOverlappy)
1. [qualityBeds.R](#qualityBedsR)
1. [sigTargets.py](#sigTargetspy)
1. [targetEnrichment.R](#targetEnrichmentR)
1. [Activity Log](#Activity%20Log)


## [Overview](https://miro.com/app/board/o9J_km5-viQ=/)
![Workflow](workflow.png)


## [betaBatch.py](betaBatch.py)
**Dependencies:** Python 2.7 and [BETA](http://cistrome.org/BETA/)

Runs BETA minus (hg38) on input bed files to generate lists of potential targets ranked by regulatory potential as it is defined in [PMC4135175](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4135175/).

## [cobinding.R](cobinding.R)
_**Note:** Designed to run as a job array with 2 cores for each job (~4GB RAM/core)_

**Dependencies:** R 4.0+, tidyverse, pheatmap, dendsort, RColorBrewer, fastcluster

Creates separate dataframes for each highly enriched or depleted TF by concatenating all the [rank normalized](#normalizepy), [peak overlapped](#peakOverlappy) bed files for that individual TF. After performing hierarchical clustering on both the rows and columns (Euclidean distance, Ward linkage) using the incredibly memory efficient [fastcluster](http://danifold.net/fastcluster.html) library, it creates cobinding heatmaps as seen in [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/). 

For each heatmap, the titular TF is the so-called "focus factor" for that graphic. The value at a given row (potential cobinding factor) and column (focus factor peak) is the rank normalized score of the potential cobinding factor peak that overlaps with the focus factor peak. If no overlapping peak exists, a score of 0 is assigned. If multiple overlapping peaks exist, the highest score among them is assigned.

By default, heatmaps are created for each TF in the top 20 by log2 fold change DE for which there is data. The rows of the heatmaps are all DE TF's that change in the same direction as the focus factor (all enriched or all depleted with quiescence).

## [deTargets.R](deTargets.R)
**Dependencies:** R 4.0+ and tidyverse

Inner joins each [combined targets file](#sigTargetspy) with differential expression data generated using DESeq2.

## [fibroblast.py](fibroblast.py)
**Dependencies:** Python 3.8+, tqdm, and pandas

Isolates all bed files for human TFs in dermal fibroblasts (though that can be configured to any other cell/tissue type) that meet the same quality criteria detailed in [qualityBeds.R](#qualitybedsR). Merges all bed files for a given TF into a single file by combining overlapping peaks into singular entries (final column shows the number of peaks combined to make a given entry).

## [greatBatch.R](greatbatch.R)
**Dependencies:** Python 3.8+, bedtools, tqdm, and pandas

Isolates all bed files for human TFs in dermal fibroblasts (though that can be configured to any other cell/tissue type) that meet the same quality criteria detailed in [qualityBeds.R](#qualitybedsR). Merges all bed files for a given TF into a single file by combining overlapping peaks into singular entries (final column shows the number of peaks combined to make a given entry).

## [greatBatch.R](greatBatch.R)
**Dependencies:** R 4.0+, rGREAT, GenomicRanges

Runs each DE TF bed file through GREAT with peaks called on unassigned scaffolds removed (chromosome names GL* or KI* as described [here](https://github.com/dpryan79/ChromosomeMappings/blob/master/GRCh38_ensembl2UCSC.txt)). Outputs separate tables of enrichment information for Biological Processes, Cellular Component, and Molecular Function. Also creates a file with several key summary graphics.

## [H4K20me3.py](H4K20me3.py)
**Dependencies:** Python 3.8+, bedtools, tqdm, and pandas

Isolates all bed files in [Cistrome database](http://cistrome.org/db/#/) corresponding to H4K20me3-targeting experiments and merges all peaks within a specified distance to create a single file of with, hopefully, all unique binding sites. The columns for the resulting file are chromosome, start, end, and number of peaks merged to create that final listing.

Then, uses [UROPA](https://www.nature.com/articles/s41598-017-02464-y#Sec2) to annotate each peak with the nearest protein coding gene. The tool also produces some nice summary graphs in a pdf file. The UROPA `json` file is configured as below:

    {
        "queries": [
            {
            "feature": "gene",
            "feature.anchor": "start",
            "distance": 3000,
            "filter.attribute": "gene_type",
            "attribute.value": "protein_coding"
            }
        ],
        "show_attributes": "gene_name",
        "gtf": "Path to gencode.v35.annotation.gtf",
        "bed": "Path to merged peak file",
        "threads": 8
    }

Additionally, this script reads the [deTargets.R](#detargetsR) output file for each DE TF and performs an inner join between those and the UROPA-annotated peak file (specifically, the one with `finalhits` in its name). It tacks on an extra column for each of these inner-joined dataframes denoting the DE TF whose targets it has used to perform this inner join. This makes it possible to distinguish which targets H4K20me3 shares with each DE TF when all of these dataframes are subsequently concatenated and written into one very long text file with the word `merged` in its name. To ease viewing and future data manipulation, the rows in this giant output file are first sorted alphabetically by DE TF name (so all targets for each DE TF are grouped together), then by location (from the first base pair of chromosome 1 to the last base pair of chromosome Y).

Lastly, the UROPA-annotated peak file is merged with the original quiescence DE data, adding log2 fold changes and DESeq2 p-adjusted values to any genes that exhibit significant (p-adj < 0.05) differential expression.

## [mergeDEandTF.R](mergeDEandTF.R)
**Dependencies:** R 4.0+ and tidyverse

Performs an inner join between a [list of human transcription factors](http://humantfs.ccbr.utoronto.ca/download.php) and differential expression data generated using DESeq2. Sorts the resulting dataframe by ascending log2 fold change and removes any rows with:  

    abs(log2FoldChange) < 1.0 OR padj > 0.05


## [normalize.py](normalize.py)
**Dependencies:** Python 3.8+, tqdm, and pandas

Reads bed files, ranks peaks, and assigns them each a normalized rank score (see section C.2.2 in [supplementary paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/bin/NIHMS541492-supplement-Supplementary_Material.pdf) to [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/)) equal to `(n-r)/(n-1)`, where `r` is the rank of each peak in the total set of `n` peaks. Output files only retain columns for the chromosome, peak start, peak end, and normalized rank.

**Note:** Output files are ranked by descending signal value (enrichment) in `$SCRATCH/tfchip/normalize/signal/` and by ascending q value (p-adjusted) in `$SCRATCH/tfchip/normalize/qval/`.


## [peakOverlap.py](peakOverlap.py)
**Dependencies:** Python 3.8+, tqdm, pandas, and bedtools

For each transcription factor, runs:

    bedtools intersect -a {factor} -b {cofactor} -loj

`factor`: the focus TF's [highest quality](#qualityBedsR), [rank normalized](#normalizepy) bed file  
`cofactor`: any highest quality, rank normalized bed file (including the focus TF's)

Each row in the output files contains the chromosome, peak start, and peak end values of a peak in `factor` and the normalized rank of the overlapping peak in `cofactor` (-1 if none found; write multiple rows if more than one found).


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

## Activity Log
1. Using [mergeDEandTF.R](#mergedeandtfr), I merged the provided differential expression data (serum starved vs proliferating cells) with a [curated list of human TFs](http://humantfs.ccbr.utoronto.ca/download.php) to create the file `SSvsP_SigTFs_080720.csv`. Only the TFs that exhibited significant differences in expression (DE TFs) were included in the final file.

2. Using [qualityBeds.R](#qualitybedsR), I isolated the bed files in the [Cistrome database](http://cistrome.org/db/#/) that corresponded to each DE TF. These files were filtered, ranked, and renamed according to specific [quality criteria](#qualitybedsR).

3. Using [betaBatch.py](#betabatchpy), I generated lists of putative targets for each high quality bed file.

4. Using [sigTargets.py](#sigtargetspy), I consolidated the BETA target predictions to create a single file for each DE TF containing all unique predictions for all high quality bed files associated with that TF. These are the files in the folder named `BETA_minus_output`.

5. Using [deTargets.R](#detargetsR), I merged the consolidated BETA target predictions for each DE TF with the differential expression data, keeping only the predicted targets that exhibited significant differential expression. These are the files in the folder named `detargets`.

6. Using [targetEnrichment.R](#targetEnrichmentR), I attempted to better characterize the enrichment behaviors of the DE TF by comparing the differential expression profiles of their targets. The resulting files are in the folder named `enrich`.

7. Using [normalize.py](#normalizepy), I assigned a rank normalized score to each peak in each high quality DE TF bed file. I generated two output directories, one for scores assigned with peaks ranked by descending signal value and another for peaks ranked by ascending p-adjusted value (see [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) for more details about how MACS2 calculates these values).

8. Using [peakOverlap.py](#peakoverlappy), I used `bedtools intersect` to search for overlapping peaks in every possible pairing of highest quality, rank normalized bed files (only one such file should exist for each DE TF). This also appended the rank normalized score of each peak to every row corresponding to a peak it overlaps.

9. Using [cobinding.R](#cobindingR), I created heatmaps for each TF that theoretically illustrate its cobinding affinity with the other DE TFs at its experimental binding loci. These heatmaps are located within the `cobinding` folder. `cobindingOLD` has these same heatmaps but as jpgs instead of pdfs. These load a lot more quickly and take up less space but have some strange visual artifacts.

10. Using [greatBatch.R](#greatbatchR), I performed GO annotation on each of the DE TF bed files using the online [GREAT](http://great.stanford.edu/public/html/) tool. These results are located in the `GO` folder.

11. Using [H4K20me3.py](#h4k20me3py), I pulled all the H4K20me3-associated bed files from the [Cistrome database](http://cistrome.org/db/#/), merged them into one file (doing my best to combine and count duplicates), and annotated them with their nearest protein coding genes (by distance to TSS). Out of curiosity, I first configured the code to only consider directly overlapping peaks as duplicates that should be merged, then tried the same pipeline with peaks within 1000bp flagged as duplicates. I then merged the list of H4K20me3-associated genes with the DE TF putative target lists to create massive text files containing information about the "target genes" this histone modification shares with the DE TFs. All of these results are located in the `H4K20me3` folder, with the direct overlap files all having a filename suffix of `0` while all the 1000bp distance outputs have a filname suffix of `1000`.

12. Using [fibroblast.py](#fibroblastpy), I pulled all the Cistrome bed files corresponding to ChIP-Seq experiments done on dermal fibroblasts. For each TF, I created a file that contains all the unique peaks called across all the skin fibroblast experiments associated with that TF. These files are located in the `skin_fibroblast` folder.

13. Using an updated version of [cobinding.R](#cobindingR), I created heatmaps for each DE TF that only include TF's whose expression changed in the same direction as the focus factor. These graphics are located in the `cobindingAll` folder. I did this again after further limiting the co-factors to those both change in the same direction and are within the top 20 by DE, creating the graphics in the `cobindingTop` folder.

14. Using an updated version of [H4K20me3.py](#h4k20me3py), I merged the UROPA-annotated H4K20me3 peak files with DE TF BETA outputs and the original DE spreadsheet to look for trends in the DE of genes located near this histone mark, whose presence is known to increase during quiescence.
