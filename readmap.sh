#!/bin/bash

# This script takes a fastq file, aligns reads using Bowtie2, and removes multi-mappers and duplicate reads.
# USAGE: sh readmap.sh <path to fastq file>

# allows module load to work
source /etc/profile.d/modules.sh

# fastq file name
fq=$1

# name outputs with basename
base=`basename $fq .fastq`
echo "Sample name is $base"

# directory with genome build
genome=${SCRATCH}/hg38.p12/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index

# make output directories
mkdir -p ${SCRATCH}/bowtie2/intermidiate_bams

# set up file names
align_out=${SCRATCH}/bowtie2/${base}_unsorted.sam
align_sorted=${SCRATCH}/bowtie2/${base}_sorted.bam
align_duped=${SCRATCH}/bowtie2/${base}_sorted_marked_dup.bam
align_filtered=${SCRATCH}/bowtie2/${base}_aln.bam

# variables of directories for cleanup
bowtie_results=${SCRATCH}/bowtie2
intermediate_bams=${SCRATCH}/bowtie2/intermidiate_bams

# load required modules
module load bowtie2
module load samtools
export PATH=$SCRATCH/gatk-4.1.9.0:$PATH

echo "Processing file $fq"

cores=$(nproc)

# Run bowtie2
bowtie2 -p ${cores} -q --local -x ${genome} -U ${fq} -S ${align_out}

# Sort BAM file by genomic coordinates
samtools sort -@ $cores -o $align_sorted $align_out

# Index sorted BAM file
samtools index -@ $cores $align_sorted

# Filter out duplicates
gatk MarkDuplicates -I $align_sorted -O $align_duped -M marked_dup_metrics.txt -QUIET True -VERBOSITY ERROR
samtools view -bF 1028 -o $align_filtered $align_duped

# Index filtered file
samtools index $align_filtered

# Move intermediate files into subdir
mv ${bowtie_results}/${base}*sorted* ${intermediate_bams}