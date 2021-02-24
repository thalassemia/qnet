# Read alignment with on very sensitive mode with soft-clipping, up to 10 multi mapping alignments
bowtie2 --very-sensitive-local -k 10 -x $SCRATCH/hg38.p12/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -p 36 -1 SS-12-3_S14_R1_001.fastq.gz -2 SS-12-3_S14_R2_001.fastq.gz | samtools view -@ 36 -u - | samtools sort -@ 36 -o ss123.bam -
samtools index -b -@ 36 ss123.bam
# Remove reads aligning to mitochondrial chromosome
samtools view -@ 36 -h ss123.bam | grep -v chrM | samtools sort -@ 36 -O bam -o ss123.rmchrm.bam
gatk --java-options "-XX:ParallelGCThreads=36 -Djava.io.tmpdir=/tmp" MarkDuplicates QUIET=true INPUT=ss123.rmchrm.bam OUTPUT=ss123.marked.bam METRICS_FILE=ss123.marked.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp
# Remove alignments with MAPQ < 30 (multi mapping alignments)
samtools view -@ 36 -h -q 30 -o ss123.rmMulti.bam ss123.marked.bam 
# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
# Retain properly paired reads -f 2
samtools view -@ 36 -h -b -F 1804 -f 2 -o ss123.filtered.bam ss123.rmMulti.bam
samtools index -b -@ 36 ss123.filtered.bam
# use --ATACshift to offset + reads by +4 bp and - reads by -5bp to account for 9bp gap between Tn5 adapter insertion
alignmentSieve --numberOfProcessors 36 --ATACshift --bam ss123.filtered.bam -o ss123.tmp.bam
# the bam file needs to be sorted and indexed again
samtools sort -@ 36 -O bam -o ss123.shifted.bam ss123.tmp.bam
samtools index -@ 36 ss123.shifted.bam
rm ss123.tmp.bam

macs3 callpeak -f BAMPE -g hs --keep-dup all -n ss123 -t ss123.shifted.bam