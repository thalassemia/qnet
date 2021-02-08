# **Box Folder Organization**

## **Table of Contents**

- [ATAC-seq](#atac-seq)
    - [Annotations](#annotations)
    - [DA-DE Correlation](#da-de-correlation)
    - [Data](#ATAC-seq-data)
- [RNA-seq](#rna-seq)
    - [Predicted Regulators](#predicted-regulators)
    - [Compare Conditions](#compare-conditions)
    - [Data](#rna-seq-data)
- [ChIP-seq](#chip-seq)
    - [Transcription Factor](#transcription-factor)
        - [Co-association](#co-association)
        - [BETA](#beta)
        - [GREAT](#great)
    - [Histone Modification](#histone-modification)
        - [Co-association](#co-association)
        - [H4K20me3](#h4k20me3)
        - [CRC](#crc)
- [RSGDREAM 2020](#rsgdream-2020)
- [Paper Notes](#paper-notes)
- [Progress Reports](#progress-reports)
- [To-Do](#to-do)


## **ATAC-seq**
### **Annotations**
Contains files for SS vs P DA regions annotated using the databases and tools described [here](atac/info.md). Also contains bar charts showing the relative frequencies of each annotation for regions with increased accessibility and those with decreased accessibility. Also contains hg38 reference sequences for DA regions.

### **DA-DE Correlation**
Contains scatter plots and inner joins representing a naive attempt to find some sort of correlation between differential expression and differential accessibility for SS vs P.

### **ATAC-seq Data**
Contains DiffBind and raw ATAC-seq output for SS vs P.


## **RNA-seq**
### **Predicted Regulators**
Contains [LISA](https://pubmed.ncbi.nlm.nih.gov/32033573/) and [ChEA3](https://pubmed.ncbi.nlm.nih.gov/31114921/) output for top 500 SS vs P DE genes (by ascending p-adjusted).
### **Compare Conditions**
Contains all output of [Comparing DE Between Quiescent Conditions](README.md#comparing-DE-between-quiescent-conditions), including CSV files of SQL-style outer joins between DE gene lists and a matrix of counts for overlaps between lists.
### **RNA-seq Data**
Contains DEseq2 output for SS vs P, CI vs P, SS vs SSR, and CI vs CIR.


## **ChIP-seq**
### **Transcription Factor**
#### **Co-association**
Contains [co-binding maps](tfchip/info.md#Visualize-Co-binding-as-Heatmap) created using various sets of SS vs P DE TFs and DE genes (e.g. upregulated, downregulated, all). Contains [interaction heatmaps](tfchip/info.md#train/interpret-ml-models) corresponding to each of these co-binding maps.
#### **BETA**
Contains [BETA output](beta/info.md#betaBatchpy) (basic and minus mode) for SS vs P DE TF ChIP-seq data from 2018 Cistrome batch download. Includes two lists of merged (deduplicated) predicted targets for each DE TF, one for all targets and a second specifically for DE targets. Also contains [heatmaps](beta/info.md#targetEnrichmentR) showing log 2 fold changes of DE targets for DE TFs at various regulatory potential cutoffs.
#### **GREAT**
Contains [GREAT output](README.md#great-go-of-chip-seq-peak-files) for all SS vs P DE TF ChIP-seq files in 2018 Cistrome batch download.

### **Histone Modification**
#### **Co-association**
Contains [co-binding maps](hmchip/info.md#create-co-binding-map) created using various sets of SS vs P DE TFs and either DE genes or DA regions.
#### [**H4K20me3**](#h4k20me3-patterns)
Contains UROPA-annotated output for merged H4K20me3 ChIP-seq data in 2018 Cistrome Batch download. Merging was performed by reading all the ChIP-seq peak files into memory, concatenating them into a single master peak list, and combining all peaks within a user-definable distance of one another (in this case, 0 and 1000 bp) using `bedtools merge`. Annotation was performed using GENCODE V29.

Also contains overlap between merged H4K20me3 annotated genes and merged SS vs P DE TF BETA-predicted target genes. Lastly, contains overlap between merged H4K20me3 annotated genes and SS vs P DE genes.
#### **CRC**  
Contains the [aligned reads](README.md#read-alignment), [ROSE output](CRC/info.md#find-putative-super-enhancers), and [CRCmapper output](CRC/info.md#create-crcs-using-crcmapper) for H3K27ac ChIP-seq data from a [cell cycle study](https://pubmed.ncbi.nlm.nih.gov/28289232/) and an [ENCODE experiment](https://www.encodeproject.org/experiments/ENCSR000APN/) on a line of adult human dermal fibroblasts. The raw sequencing reads for G0/G1, M, and S phase H3K27ac ChIP-seq were downloaded using [SRA Explorer](https://sra-explorer.info/) and the search term SRP098814. Peaks were subsequently called using MACS3. By contrast, prealigned reads (bowtie2, hg38) from isogenic replicate 2 (ENCFF754VWN) were downloaded for the ENCODE cell line. Similarly, pre-called replicated peaks (ENCFF398SWW) were downloaded for this experiment.


## **RSGDREAM 2020**
Contains the abstract, poster, and recorded five-minute talk submitted to the virtual 2020 ISCB Regulatory and Systems Genomics with DREAM Challenges Conference.


## **Paper Notes**
Contains notes on relevant papers for future reference.

## **Progress Reports**
Contains lab meeting presentations and other progress reports.

## **To-Do**
Contains the most up-to-date lists of future plans and extensions.