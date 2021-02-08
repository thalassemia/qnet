# ATAC-seq Annotation

## Table of Contents
1. [Annotate Differentially Accessible (DA) Regions](#annotate-da-regions)
1. [Annotate DA Regions with Candidate Cis-Regulatory Elements (cCRE)](#compare-with-ccres)
1. [Create Bar Charts for Annotation Frequencies](#create-bar-charts-for-annotation-frequencies)
1. [Fetch Reference Sequences for DA Regions](#fetch-reference-sequences)
1. [To-Do](#to-do)

---

## [Annotate DA Regions](anno.r)

Uses ChIPseeker and GENCODE V29 to annotate ATAC-seq DA regions. Filters for DA regions within 3 kb of a promoter and plots DA log 2 fold change against RNA-seq DE log 2 fold change. Only DA and DE regions with p-adjusted < 0.05 are considered.

## [Compare with cCREs](overlaps.R)

Compares ChIPseeker output from above to three cCRE annotation databases:
1. [HoneyBadger2](https://personal.broadinstitute.org/meuleman/reg2map/HoneyBadger2_release/): DNaseI hypersensitivity regions across 53 epigenomes with -log<sub>10</sub>(p)>2 (***hg19***)
1. [Roadmap 15 Mark ChromHMM](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html): E126 (adult human dermal fibroblast), hg38 liftover
1. [ENCODE SCREEN](https://screen.encodeproject.org/): integrates DNase, H3K4me3, H3K27ac, and CTCF data to identify cCREs

Specifically, runs the following command:

    bedtools intersect -wao -a <chipseeker> -b <ccre>

## [Create Bar Charts for Annotation Frequencies](charts.py)

Runs `bedtools merge` to collapse all cCRE annotations for a given DA region into a single line:

    chr1 3454156 3454656 ... pELS,CTCF-bound
    chr1 3454156 3454656 ... PLS,CTCF-bound

    bedtools merge ...

    chr1 938315 938815 ... pELS,CTCF-bound,PLS,CTCF-bound

Remove duplicate annotations (e.g. CTCF-bound above) and separate DA regions into two lists: one for increased accessibility and another for decreased accessibility with quiescence. Calculate frequency of each annotation combination and plot annotation frequencies as horizontal bar chart (one bar for up and another for down).

**Special Notes**
1. ChIPseeker: introns and exons by default have very long tags (e.g. `Intron (ENST00000370056.8/ENSG00000134215.15, intron 20 of 26)`). These were all collapsed into two discrete categories: `Intron` and `Exon`.
1. ChromHMM: annotations are somewhat esoteric (e.g. `5_TxWk`). These were replaced by more parsable labels (e.g. Weak transcription)

## [Fetch Reference Sequences](getSeq.r)

Uses `BSgenome.Hsapiens.UCSC.hg38` from Bioconductor to fetch list of reference sequences for each DA region.