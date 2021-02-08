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