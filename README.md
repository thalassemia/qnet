# **Genomic Exploration of Quiescence**
_**Using a mix of in-house and public datasets, we aim to extrapolate meaningful information about the regulatory circuitry underlying the maintenance of cellular quiescence.**_

## **Table of Contents**

**TF Co-association Pipeline**
1. [Overview](tfchip/Poster.pdf)
1. [Filter ChIP-Seq Data](tfchip/info.md#filter-chip-seq-data)
1. [Preprocess Data](tfchip/info.md#preprocess-data)
1. [Train/Interpret ML Models](tfchip/info.md#train/interpret-ml-models)
1. [Analyze Results](tfchip/analyze/info.md)

---

**ATAC-seq Annotation**
1. [Annotate Differentially Accessible (DA) Regions](atac/info.md#annotate-da-regions)
1. [Annotate DA Regions with Candidate Cis-Regulatory Elements (cCRE)](atac/info.md#Compare-with-cCRE-databases)
1. [Create Bar Charts for Annotation Frequencies](atac/info.md#create-bar-charts-for-annotations)
1. [Fetch DNA Sequences for DA Regions](atac/info.md#fetch-reference-sequences)

---

**Histone Modifications**
1. [Filter ChIP-seq Data](hmchip/info.md#filter-chip-seq-data)
1. [Preprocess Data](hmchip/info.md#preprocess-data)
1. [Create Co-binding Map](hmchip/info.md#create-co-binding-map)

---

**Core Regulatory Circuits (CRCs)**
1. [Convert Ensembl to HGNC IDs](CRC/info.md#convert-ensembl-to-hgnc-ids)
1. [Find Putative Super Enhancers](CRC/info.md#find-putative-super-enhancers)
1. [Create CRCs using CRCmapper](CRC/info.md#create-crcs-using-crcmapper)

---

**Miscellaneous**
1. [Determine DE TF Targets](beta/info.md)
1. [GREAT GO of ChIP-Seq Peak Files](#great-go-of-chip-seq-peak-files)
1. [H4K20me3 Patterns](#h4k20me3-patterns)
1. [Comparing DE Between Quiescent Conditions](#comparing-DE-between-quiescent-conditions)
1. [Read Alignment](#read-alignment)
1. [Box Folder Organization](box.md)

## [GREAT GO of ChIP-Seq Peak Files](greatBatch.R)

Runs each DE TF bed file through GREAT with peaks called on unassigned scaffolds removed (chromosome names GL* or KI* as described [here](https://github.com/dpryan79/ChromosomeMappings/blob/master/GRCh38_ensembl2UCSC.txt)). Outputs separate tables of enrichment information for Biological Processes, Cellular Component, and Molecular Function. Also creates a file with several key summary graphics.

## [H4K20me3 Patterns](H4K20me3.py)

Isolates all bed files in [Cistrome database](http://cistrome.org/db/#/) corresponding to H4K20me3-targeting experiments and merges all peaks within a specified distance to create a single file of with, hopefully, all unique H4K27me3 deposition sites. The columns for the resulting file are chromosome, start, end, and number of peaks merged to create each final entry.

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
        "gtf": "Path to gencode.v29.gtf",
        "bed": "Path to merged peak file",
        "threads": 36
    }

Additionally, this script reads the [deTargets.R](#detargetsR) output file for each DE TF and performs an inner join between those and the UROPA-annotated peak file (specifically, the one with `finalhits` in its name). It tacks on an extra column for each of these inner-joined dataframes denoting the DE TF whose targets it has used to perform this inner join. This makes it possible to distinguish which targets H4K20me3 shares with each DE TF when all of these dataframes are subsequently concatenated and written into one very long text file with the word `merged` in its name. To ease viewing and future data manipulation, the rows in this giant output file are first sorted alphabetically by DE TF name (so all targets for each DE TF are grouped together), then by location (from the first base pair of chromosome 1 to the last base pair of chromosome Y).

Lastly, the UROPA-annotated peak file is merged with the original quiescence DE data, adding log2 fold changes and DESeq2 p-adjusted values to any genes that exhibit significant (p-adj < 0.05) differential expression.

## [Read Alignment](readmap.sh)

**Note**: Requires gatk, samtools, and bowtie2 (will replace with BWA)

Script to automate filtering and alignment of raw Illumina sequencing reads. Uses GRCh38.p12 (same patch version as GENCODE V29) bowtie index from [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/GRCh38_major_release_seqs_for_alignment_pipelines/).

## [Comparing DE Between Quiescent Conditions](newCond.py)
Filters all DE gene lists (DEseq2) for p-adjusted <= 0.05 and abs(log2FoldChange) >= 1. Creates secondary lists for each primary list that only contain DE TFs. Performs SQL-style outer joins between each possible pairing of DE lists (does NOT mix-and-match between gene and TF lists). Also includes two matrices of counts (one for all DE genes and another for DE TFs) with the following format:

<center>

&nbsp;       | SS vs P | CI vs P | SS vs SSR | CI vs CIR
-------------|---------|---------|-----------|----------
**SS vs P**  |Total    |Same Dir.|Same Dir.  |Same Dir.
**CI vs P**  |Diff Dir.|Total    |Same Dir.  |Same Dir.
**SS vs SSR**|Diff Dir.|Diff Dir.|Total      |Same Dir.
**CI vs CIR**|Diff Dir.|Diff Dir.|Diff Dir.  |Total

</center>

&nbsp;  
**SS** = serum starved, **CI** = contact inhibited, **P** = proliferating, **SSR** = serum starved restimulated, **CIR** = contact inhibited restimulated, **Total** = count for that DE list, **Same Dir.** = count of overlapping entries with the same fold change sign in both commpared lists, **Diff Dir.** = count of overlapping entries with different fold change signs in both compared lists

## Activity Log
1. Using [mergeDEandTF.R](getBeds/mergeDEandTF.R), I merged the provided differential expression data (serum starved vs proliferating cells) with a [curated list of human TFs](http://humantfs.ccbr.utoronto.ca/download.php) to create the file `SSvsP_SigTFs_080720.csv`. Only the TFs that exhibited significant differences in expression (DE TFs) were included in the final file.

2. Using [qualityBeds.R](getBeds/qualityBeds.R), I isolated the bed files in the [Cistrome database](http://cistrome.org/db/#/) that corresponded to each DE TF. These files were filtered, ranked, and renamed according to specific [quality criteria](#qualitybedsR).

3. Using [betaBatch.py](beta/betaBatch.py), I generated lists of putative targets for each high quality bed file.

4. Using [sigTargets.py](beta/sigTargets.py), I consolidated the BETA target predictions to create a single file for each DE TF containing all unique predictions for all high quality bed files associated with that TF. These are the files in the folder named `BETA_minus_output`.

5. Using [deTargets.R](beta/deTargets.R), I merged the consolidated BETA target predictions for each DE TF with the differential expression data, keeping only the predicted targets that exhibited significant differential expression. These are the files in the folder named `detargets`.

6. Using [targetEnrichment.R](beta/targetEnrichment.R), I attempted to better characterize the enrichment behaviors of the DE TF by comparing the differential expression profiles of their targets. The resulting files are in the folder named `enrich`.

7. Using normalize.py (consolidated into files labelled `prep`), I assigned a rank normalized score to each peak in each high quality DE TF bed file. I generated two output directories, one for scores assigned with peaks ranked by descending signal value and another for peaks ranked by ascending p-adjusted value (see [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) for more details about how MACS2 calculates these values).

8. Using [peakOverlap.py](preprocess/peakoverlap.py), I used `bedtools intersect` to search for overlapping peaks in every possible pairing of highest quality, rank normalized bed files (only one such file should exist for each DE TF). This also appended the rank normalized score of each peak to every row corresponding to a peak it overlaps.

9. Using [cobindO.R](cobindHeatmaps/cobindO.R), I created heatmaps for each TF that theoretically illustrate its cobinding affinity with the other DE TFs at its experimental binding loci. These heatmaps are located within the `cobinding` folder. `cobindingOLD` has these same heatmaps but as jpgs instead of pdfs. These load a lot more quickly and take up less space but have some strange visual artifacts.

10. Using [greatBatch.R](#greatbatchR), I performed GO annotation on each of the DE TF bed files using the online [GREAT](http://great.stanford.edu/public/html/) tool. These results are located in the `GO` folder.

11. Using [H4K20me3.py](#h4k20me3py), I pulled all the H4K20me3-associated bed files from the [Cistrome database](http://cistrome.org/db/#/), merged them into one file (doing my best to combine and count duplicates), and annotated them with their nearest protein coding genes (by distance to TSS). Out of curiosity, I first configured the code to only consider directly overlapping peaks as duplicates that should be merged, then tried the same pipeline with peaks within 1000bp flagged as duplicates. I then merged the list of H4K20me3-associated genes with the DE TF putative target lists to create massive text files containing information about the "target genes" this histone modification shares with the DE TFs. All of these results are located in the `H4K20me3` folder, with the direct overlap files all having a filename suffix of `0` while all the 1000bp distance outputs have a filname suffix of `1000`.

12. Using an updated version of [cobindO.R](cobindHeatmaps/cobindO.R), I created heatmaps for each DE TF that only include TF's whose expression changed in the same direction as the focus factor. These graphics are located in the `cobindingAll` folder. I did this again after further limiting the co-factors to those both change in the same direction and are within the top 20 by DE, creating the graphics in the `cobindingTop` folder.

13. Using an updated version of [H4K20me3.py](#h4k20me3py), I merged the UROPA-annotated H4K20me3 peak files with DE TF BETA outputs and the original DE spreadsheet to look for trends in the DE of genes located near this histone mark, whose presence is known to increase during quiescence.

14. I used the top 500 DE genes (by ascending p-adj.) as inputs to ChEA3 and Lisa. Using [chea_lisa_beds.R](getBeds/chea_lisa_beds.R), I found all Cistrome peak files corresponding to the top 200 TFs for the output of each webtool. Using [peakAnno.R](preproces/peakAnno.R), I annotated these peak files using ChIPseeker. Using [chea_lisa_prep.py](preprocess/chea_lisa_prep.py), I preprocessed these annotated peak files as described [here](preprocess/info.md).

15. Using [peakAnno.R](preproces/peakAnno.R), I annotated the peak files for the DE TFs (from step 2) using ChIPseeker. Using [up_prep.py](preprocsess/up_prep.py), I preprocessed these peak files as described [here](preprocess/info.md).

16. Using the files in the [ml folder](ml/info.md), I trained and saved 100 LightGBM models and their corresponding SHAP outputs for each of my three input TF lists (ChEA3, Lisa, and upregulated TFs).

17. Using [acc_stats.py](analyze/acc_stats.py), I determined the average training and test accuracies of models trained on my three input TF lists. 

18. Using the files in the [interactHeatmaps](interactHeatmaps/info.md) and [cobindHeatmaps](cobindHeatmaps/info.md) folders, I created heatmap visualizations of co-binding matrices and SHAP interaction value matrices.