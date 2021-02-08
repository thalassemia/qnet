library(GenomicFeatures)
library(regioneR)
library(rtracklayer)
library(Rsamtools)
library(data.table)
library(ggplot2)

#my attempt at implementing rose algorithm

#wrapper for bedtools coverage function
bedTools.cov<-function(bed, bam)
{
  #create temp files
  a.file = tempfile()
  out = tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write bed formatted dataframes to tempfile
  write.table(bed,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)

  #create the command string and call the command using system()
  command = paste("bedtools coverage -a", a.file, "-b", bam, "-counts >", out)
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F)
  unlink(a.file);
  unlink(out)
  return(res)
}

#wrapper for bedtools merge function (within 12.5kb)
bedTools.merge<-function(bed)
{
  #create temp files
  a.file = tempfile()
  out = tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write bed formatted dataframes to tempfile
  write.table(bed,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)

  #create the command string and call the command using system()
  command = paste("bedtools merge -i", a.file, "-d 12500 >", out)
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F)
  unlink(a.file);
  unlink(out)
  return(res)
}

#create and save Gencode V29 TxDb annotation database
#gencode <- system("echo $SCRATCH/data/gencode.v29.annotation.gtf.gz", intern = T)
#gencode <- makeTxDbFromGFF(gencode, format = "gtf")
#saveDb(gencode, system("echo $SCRATCH/data/gencode_v29.sqlite", intern = T))

gencode <- loadDb(system("echo $SCRATCH/data/gencode_v29.sqlite", intern = T))

#find h3k27ac peaks that do not fall within +/-2000bp of a TSS
pr <- promoters(gencode, upstream = 2000, downstream = 2000)
peaks <- fread(system("echo $SCRATCH/cellcyc/macs3/s_peaks.broadPeak", intern = T))
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(peaks$V2, end = peaks$V3, names = peaks$V4))
enhancers <- peaks[!overlapsAny(peaks, pr)]
enhancers <- sort(enhancers)
df <- data.frame(seqnames = seqnames(enhancers), names = names(enhancers), 
                from = c(rep(".", length(enhancers))), 
                starts=start(enhancers), ends = end(enhancers), 
                score = c(rep(".", length(enhancers))), strands = c(rep(".", length(enhancers))), 
                rand = c(rep(".", length(enhancers))), id = names(enhancers))
write.table(df, file=system("echo $SCRATCH/cellcyc/enhancers/m_enhancers.gff", intern = T), quote=F, sep="\t", row.names=F, col.names=F)
#merge enhancers within 12.5kb into enhancer domains
edomains <- bedTools.merge(bed = enhancers)

#extract experimental and control read coverages at enhancer domains, normalize to CPM, and subtract
experimental <- bedTools.cov(bed = edomains, bam = system("echo $SCRATCH/bowtie2/SRR5227994_GSM2476263_MCF7_H3K27ac_M_phase_Homo_sapiens_ChIP-Seq_aln.bam", intern = T))
#control <- bedTools.cov(bed = edomains, bam = system("echo $SCRATCH/h3k27ac/control.bam", intern = T))
normed <- experimental
#normed[,4] <- (experimental[,4]/25324056 - control[,4]/41347323) * 1E6
normed <- normed[order(normed[,4]),]
normed$rank <- 1:nrow(normed)

#plot graph of rank to enrichment
pdf(system("echo $SCRATCH/cellcyc/enhancers/m_plot.pdf", intern = T))
ggplot(normed, aes(x = rank, y = V4)) + 
  geom_point()
dev.off()