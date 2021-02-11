# **TF Co-association Pipeline**

## **Table of Contents**
### Filter ChIP-seq Data
1. [Get DE TF List](#get-de-tf-list)
1. [ENCODE ChIP-seq Data](#encode-chip-seq-data)
1. [Cistrome ChIP-seq Data](#de-tfs-chip-seq-data)

### Preprocess Data
1. [Annotate Peaks](#annotate-peaks)
1. [DE Gene Prep](#de-gene-prep)
1. [OLD: TF-centric Prep (Not updated)](#tf-centric-approach)
1. [Visualize Co-binding as Heatmap](#Visualize-Co-binding-as-Heatmap)
1. [OLD: TF-centric Co-binding (Not updated)](#old:-tf-centric-co-binding)

### Train/Interpret ML Models
1. [LightGBM and SHAP](#lightgbm-and-shap)
1. [Interaction Heatmaps](#interaction-heatmaps)


## **Filter ChIP-seq Data**
### [**Get DE TF List**](getBeds/mergeDEandTF.R)

Uses a [list of human transcription factors](http://humantfs.ccbr.utoronto.ca/download.php) and DESeq2 output to identify differentially expressed TFs by Ensembl gene ID. Sorts the resulting dataframe by ascending log2 fold change, keeping only rows with:  

    abs(log2FoldChange)>=1.0 AND padj>=0.05

### [**ENCODE ChIP-seq Data**](getBeds/encode.py)

Fetches IDR thresholded peak files (preferably "optimal", "psuedoreplicated", or "conservative"; see [Outputs](https://www.encodeproject.org/chip-seq/transcription_factor/#outputs)) for each [ENCODE](https://www.encodeproject.org/) experiment targeting a DE TF. Filters out experiments with extremely low read depth or special treatment. Compiles table of information about each downloaded file, including experiment ID, target TF, cell line, etc.

### [**Cistrome ChIP-seq Data**](getBeds/qualityBeds.R)

For each transcription factor that exhibits [significant differential expression](#get-de-tf-list), this reads the [Cistrome ChIP-Seq](http://cistrome.org/db) index file to identify the relevant bed files. Then, it uses a mix of [ENCODE](https://www.encodeproject.org/data-standards/terms/) and [Cistrome](http://cistrome.org/db/#/about) ChIP-Seq quality guidelines to filter and rank all bed files. 

Specifically, bed files must at least meet the following criteria:

    UniquelyMappedRatio>=0.5 AND fastQC>=25 AND PBC>=0.5 AND PeaksUnionDHSRatio>=0.7 

The passing files are then sorted (first by descending `FRiP` then by descending `PeaksFoldChangeAbove10`) before being copied to a new directory with the following structure:

    $SCRATCH/cobinding/beds/{TF}/{rank}.bed

Provides list of DE TFs for which neither Cistrome nor ENCODE had ChIP-seq data. Appends information about copied Cistrome files to metadata for ENCODE files.


## **Preprocess Data**

### [**Annotate Peaks**](preprocess/peakAnno.R)
For each bed file in a given directory, assigns each peak to its nearest transcript annotation using ChIPseeker (promoter set to TSS &#xb1; 3000bp) and GENCODE V29.

### [**DE Gene Prep**](preprocess/prep.py)
Performs five main steps in order as listed below:

1. Reads annotated bed files, ranks peaks by **ascending** signal/significance, and assigns them each a normalized rank score (see section C.2.2 in [supplementary paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/bin/NIHMS541492-supplement-Supplementary_Material.pdf) to [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/)) equal to `r/n`, where `r` is the rank of each peak in the total set of `n` peaks.

1. Reads normalized bed files and filters for peaks annotated to promoters of significant DE genes (matching is performed on stable Ensembl gene IDs). Keeps the highest scoring peak in the case of multiple peaks being annotated to a single DE gene ID. DE gene promoters with no corresponding ChIP-seq peak are automatically assigned a score of 0.

1. Retrieves the ranks of the Cistrome-derived bed files (see [Cistrome ChIP-seq Data](#cistrome-ChIP-seq-Data)), simultaneously assigning ENCODE-derived bed files with a rank of 0 to give them preference over all others when available.

1. Using the highest-ranked file for each TF, compiles normalized binding scores into a **co-binding map** with the following format, where S<sub>n,m</sub> represents the binding score of the m<sup>th</sup> TF at the n<sup>th</sup> DE gene promoter:
<center>

&nbsp;     | TF 1            | TF 2            | ...
-----------|-----------------|-----------------|----
Promoter 1 | S<sub>1,1</sub> | S<sub>1,2</sub> | ...
Promoter 2 | S<sub>2,1</sub> | S<sub>2,2</sub> | ...
...        | ...             | ...             | ...

</center>
&nbsp;  

5. Repeats step 4 using second-ranked file for each TF (reuse first if there is none) to create a test set co-binding map.

### [**OLD**: TF-centric Prep](preprocess/peakOverlap.py)

For each transcription factor, runs:

    bedtools intersect -a {factor} -b {cofactor} -loj

`factor`: the focus TF's [highest quality](#qualityBedsR), [rank normalized](#normalizepy) bed file  
`cofactor`: any highest quality, rank normalized bed file (including the focus TF's)

Each row in the output files contains the chromosome, peak start, and peak end values of a peak in `factor` and the normalized rank of the overlapping peak in `cofactor` (-1 if none found; write multiple rows if more than one found).

### [**Visualize Co-binding as Heatmap**](cobindMaps/cobindMap.r)
Creates a heatmap visualization of all co-binding matrices (see [DE Gene Prep](#DE-gene-prep)) in a given directory. For each heatmap, hierarchical clustering is performed on both the rows and columns using Euclidean distance and Ward linkage. The clustering dendrograms are rearranged using [dendsort](https://cran.r-project.org/web/packages/dendsort/index.html) to get the final output.

**Note**: For reasons unknown, dendsort fails with the error `C stack usage 79xxxxx is too close to the limit` when running this code on `down_train.csv`. Thus, the co-binding heatmap for this one file does NOT have its rows/columns arranged using dendsort.

### [**OLD**: TF-centric Co-binding](cobindMaps/cobindOld.R)
Creates separate dataframes for each highly enriched or depleted TF by concatenating all the rank-normalized, peak-overlapped bed files for that individual TF. After performing hierarchical clustering on both the rows and columns (Euclidean distance, Ward linkage) using the memory-efficient [fastcluster](http://danifold.net/fastcluster.html) library, it creates cobinding heatmaps as seen in [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/). 

For each heatmap, the titular TF is the so-called "focus factor" for that graphic. The value at a given row (potential cobinding factor) and column (focus factor peak) is the rank normalized score of the potential cobinding factor peak that overlaps with the focus factor peak. If no overlapping peak exists, a score of 0 is assigned. If multiple overlapping peaks exist, the highest score among them is assigned.

By default, heatmaps are created for each TF in the top 20 by log2 fold change DE for which there is data. The rows of the heatmaps are all DE TF's that change in the same direction as the focus factor (either all enriched or all depleted with quiescence).


## **Train/Interpret ML Models**

### [**LightGBM and SHAP**](ml.py)
Takes two co-binding maps as input: a training set and a test set. For the training set, treats each row (e.g. binding scores at each DE gene promoter) as a sample and each column (e.g. DE TF) as a feature. Assigns all samples in the input with a label of `1`. Independently shuffles each column to throughly disrupt any dependencies between features, creating a second set of samples that are all assigned a label of `0`.

Trains 100 LightGBM classification models on these training set-derived samples and performs validation with the input test set using the following parameters:

    'objective': 'binary'
    'metric': 'binary_logloss'
    early_stopping_rounds = 50

Calculates train and test set accuracies for each model (saved in Hoffman2 job output file).

Saves the resulting LightGBM model, the training set, the test set, the SHAP interaction values, the SHAP values, the list of features, and the SHAP explainer object in a python shelf for future access.

### [**Interaction Heatmaps**](interactMaps/info.md)

Creates two heatmaps of SHAP interaction values. To create the heatmap of averages, the SHAP interaction values for each TF pair are first averaged across all positive samples for a given model, then these averages were averaged again over all models. To create the heatmap of maximums, the SHAP interaction values at each positive sample are averaged across all models and the maximum (across all samples) is taken for each TF pair.

**Example:** For a model with 200 TFs (features) and 1000 DE genes (2000 samples after accounting for negative set), we have `200 * 200 * 2000` SHAP interaction values organized into 200-by-200 matrices for each of the 2000 samples. If we have 100 such models, we can consolidate all of these SHAP interaction values into a 4D array with the dimensions `(100, 2000, 200, 200)`. 

- All negative samples are discarded, leaving an array with the dimensions `(100, 1000, 200, 200)`.
- Heatmap of averages: Take the average across all samples to get an array with the dimensions `(100, 200, 200)`. Take the average across all models to get an array with the dimensions `(200, 200)`.
- Heatmap of maximums: Take the average across all models to get an array with the dimensions `(1000, 200, 200)`. Take the maximum across all samples to get an array with the dimensions `(200, 200)`.