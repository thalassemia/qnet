# **Box Folder Organization**

## **Table of Contents**

- [ATAC-seq](#atac-seq)
    - [Annotations](#annotations)
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
- [Meeting Notes](#meeting-notes)
- [To-Do](#to-do)

---

## **ATAC-seq**
### **Annotations**
Contains files for SS vs P DA regions annotated using the databases and tools described [here](atac/info.md). Also contains bar charts showing the relative frequencies of each annotation for regions with increased accessibility and those with decreased accessibility.

### **ATAC-seq Data**
Contains DiffBind and raw ATAC-seq output for SS vs P. Also contains hg38 reference sequences for DA regions.

---

## **RNA-seq**
### **Predicted Regulators**
Contains [LISA](https://pubmed.ncbi.nlm.nih.gov/32033573/) and [ChEA3](https://pubmed.ncbi.nlm.nih.gov/31114921/) output for top 500 SS vs P DE genes (by ascending p-adjusted).
### **Compare Conditions**
Contains all output of [newCond.py](README.md#comparing-DE-between-quiescent-conditions), including csv files of SQL-style outer joins between each possible pairing of quiescence DE gene and TF lists. Also includes two matrices of counts (one for all DE genes and another specifically for DE TFs) with the following format:

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

### **RNA-seq Data**
Contains DEseq2 output for SS vs P, CI vs P, SS vs SSR, and CI vs CIR.

---

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
Contains [co-binding maps](hmchip/info.md#create-co-binding-map) created using various sets of SS vs P DE TFs and SS vs P DE genes (e.g. upregulated, downregulated, all).
#### [**H4K20me3**](#h4k20me3-patterns)
Contains UROPA-annotated output for merged H4K20me3 ChIP-seq data in 2018 Cistrome Batch download. Merging was performed by reading all the ChIP-seq peak files into memory, concatenating them into a single master peak list, and combining all peaks within a user-definable distance of one another (in this case, 0 and 1000 bp) using `bedtools merge`. Annotation was performed using GENCODE V35.

Also contains overlap between merged H4K20me3 annotated genes and merged SS vs P DE TF BETA-predicted target genes. Lastly, contains overlap between merged H4K20me3 annotated genes and SS vs P DE genes.
#### **CRC**  
Contains the [aligned reads](README.md#read-alignment), [ROSE output](CRC/info.md#find-putative-super-enhancers), and [CRCmapper output](CRC/info.md#create-crcs-using-crcmapper) for H3K27ac ChIP-seq data from a [cell cycle study](https://pubmed.ncbi.nlm.nih.gov/28289232/) and an [ENCODE experiment](https://www.encodeproject.org/experiments/ENCSR000APN/) on a line of adult human dermal fibroblasts. The raw sequencing reads for G0/G1, M, and S phase H3K27ac ChIP-seq were downloaded using [SRA Explorer](https://sra-explorer.info/) and the search term SRP098814. By contrast, prealigned reads (BWA, hg38) from isogenic replicate 2 (ENCFF680YEV) were downloaded for the ENCODE cell line.

---

## **RSGDREAM 2020**
Contains the abstract, poster, and recorded five-minute talk submitted to the virtual 2020 ISCB Regulatory and Systems Genomics with DREAM Challenges Conference.

---

## **Paper Notes**
Contains notes on relevant papers for future reference.

---

## **Meeting Notes**
Contains notes on things discussed in weekly progress report meetings.

---

## **To-Do**
Contains the most up-to-date lists of future plans and extensions.