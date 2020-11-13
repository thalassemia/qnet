# Preprocess Data

## Table of Contents
1. [Annotate Peaks](#annotate-peaks)
1. [TF-centric Approach](#tf-centric-approach)
1. [DE Gene Approach](#de-gene-approach)

## [Annotate Peaks](peakAnno.R)
**Dependencies:** R 4.0+, tidyverse, ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, doParallel

Annotates bed files with ChIPseeker (TSS +/- 3000bp) using hg38 genome.

## [TF-centric Approach](peakOverlap.py)
**Dependencies:** Python 3.8+, tqdm, pandas, and bedtools

For each transcription factor, runs:

    bedtools intersect -a {factor} -b {cofactor} -loj

`factor`: the focus TF's [highest quality](#qualityBedsR), [rank normalized](#normalizepy) bed file  
`cofactor`: any highest quality, rank normalized bed file (including the focus TF's)

Each row in the output files contains the chromosome, peak start, and peak end values of a peak in `factor` and the normalized rank of the overlapping peak in `cofactor` (-1 if none found; write multiple rows if more than one found).

## [DE Gene Approach](up_prep.py)
**Dependencies:** Python 3.8+, tqdm, pandas 1.1.3, numpy 1.19.2, json

Reads annotated bed files, ranks peaks, and assigns them each a normalized rank score (see section C.2.2 in [supplementary paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/bin/NIHMS541492-supplement-Supplementary_Material.pdf) to [PMC4154057](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154057/)) equal to `(n-r)/(n-1)`, where `r` is the rank of each peak in the total set of `n` peaks. Output files only retain columns for the chromosome, peak start, peak end, and normalized rank.

Filters for bed files for peaks that are annotated to a significantly differentially-expressed gene (p-adj. < 0.05). Sorts bed files by descending FRiP.

Creates a co-binding matrix consolidating all the normalized peak scores for the files with the highest FRiP for each DE TF. The value at a given row (each of which represents a bin around a DE gene) and column (each of which represents a TF) is the rank normalized score of the TF peak that overlaps with the DE gene bin. If no overlapping peak exists, a score of 0 is assigned.