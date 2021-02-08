# Core Regulatory Circuits (CRCs)

## Table of Contents
1. [Convert Ensembl to HGNC IDs](#convert-ensembl-to-hgnc-ids)
1. [Find Putative Super Enhancers](#find-putative-super-enhancers)
1. [Create CRCs Using CRCmapper](#create-crcs-using-crcmapper)

---

## [Convert Ensembl to HGNC IDs](enst2hgnc.py)
Create a two column table formatted as below for use with CRCmapper.

Ensembl Transcript ID | HGNC Gene Symbol 
----------------------|-----------------
ENST00000428771.6     | HES4
ENST00000511072.5     | PRDM16
...                   | ...

## [Find Putative Super Enhancers](./rose/)
_**Note**_: Needs environment that has both R and python installed. 

With H3K27ac aligned reads (preferably both experimental and control) and MACS2 peaks as input, uses [ROSE](https://bitbucket.org/young_computation/rose/src/master/) to identify enhancers ans super enhancers (SEs) using the methodology described in the Extended Experimental Procedures section of these papers: [23582323](https://pubmed.ncbi.nlm.nih.gov/23582323/) and [23582322](https://pubmed.ncbi.nlm.nih.gov/23582322/). Briefly:

1. MACS2 peaks outside of promoter regions (defined to be TSS &#xb1; 2.5 kb) are known as _constituent enhancers_. Constituent enhancers within 12.5 kb of one another are stitched together to form _active enhancers_.
1. Aligned reads are extended by 200 bp on their 3' ends.
1. For each constituent enhancer, the number of extended reads overlapping each nucleotide position are are summed and normalized by total number of reads (rpm) and enhancer size (rpm/bp). The normalized control densities are subtracted from the normalized experimental densities to get true read densities.
1. Read densities of active enhancers are calculated by summing the true read densities of their constituent enhancers.
1. Active enhancers are ranked by ascending read density. This rank is plotted against read density. All active enhancers falling to the right of the point on the resulting graph with a slope of 1 are calssified as _super enhancers_.

## [Create CRCs Using CRCmapper](./crc/)
Uses H3K27ac aligned reads, MACS2 peaks, and a ROSE enhancer table as inputs to identify loops of self-regulating transcription factors associated with super enhancers. The exact methodology is detailed in this paper: [26843070](https://pubmed.ncbi.nlm.nih.gov/26843070/). Briefly:

1. ROSE-identified super enhancers (SEs) are assigned to their closest expressed transcripts (transcript expression estimated using H3K27ac signal by default). See [25303531](https://pubmed.ncbi.nlm.nih.gov/25303531/) for evidence that most enhancers (super or otherwise) likely regulate their most proximal genes.
1. SE-assigned TFs are identified and a motif search (FIMO, see [21330290](https://pubmed.ncbi.nlm.nih.gov/21330290/)) is conducted on the sequences of SE constituents (extended by 500 bp on each side) using the PWM of their corresponding TF. 
1. SE-assigned TFs with constituents containing 3+ instances of their motif are termed _auto-regulated TFs_ (ATFs). Motif searches are conducted on all SEs assigned to ATFs using the PWMs of all ATFs to identify fully interconnected auto-regulatory loops, which are outputted as candidate CRCs.

***Key Modifications***

1. Downloaded "comprehensive" GENCODE V29 annotation file from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) to keep analysis consistent. Required a few tweaks to the code and the replacement of the pre-bundled refseq-to-HGNC table with an [Ensembl-to-HGNC](#convert-ensembl-to-hgnc-ids) table.
1. Implemented multiprocessing to dramatically speed up the mapping of aligned reads to enhancer regions.

## **To-Do**
1. Quantify transcript expression using in-house ATAC-seq or RNA-seq data.
1. Prepare to analyze in-house CUT&Tag data for histone modifications.